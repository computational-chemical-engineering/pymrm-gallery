# Prompt for the next session

Copy everything in the block below as the first message of a fresh session.

---

```
Continue the pymrm-gallery. Repo ~/Code/pymrm_suite/pymrm-gallery, venv ~/Code/pymrm_suite/.venv.

Operate per docs/playbooks/coordinator.md: read docs/handoff.md header + git log -8,
then work ONE CASE AT A TIME — builder → verifier → fixer → integrate → commit — and
only start the next once that one is pushed. Do NOT run builders in parallel: that
lost a whole wave to a usage limit on 2026-08-07. Dispatch with the 10-line template
in docs/playbooks/README.md (the playbooks carry the rules — do not restate them) and
set model: per its tier table (verifiers never downgraded). Per-source traps are in
docs/papers-inventory.yaml; integrate with scripts/integrate_case.py; pull --rebase
before every push; never git add -A.

Suggested order — all T0/P1, sources on disk, each with published siblings to reuse
and cross-check against (reading those siblings is required, not optional):
  J1.4 IAST        — Myers & Prausnitz 1965. Completes the adsorption trio with the
                     just-published J1.1 (Langmuir) and J1.3 (BET). J1.1 found that
                     Langmuir's own Case VI IS the BET isotherm; do not restate it.
  B1.3 Bischoff    — completes the modulus quartet; B1.1 (Thiele), B1.2 (Aris) and
                     B1.5 (Weisz–Hicks) are all published and reusable.
  C1.2 Eley–Rideal — C1.1 (LHHW) is published, so import the Langmuir–Hinshelwood law
                     rather than re-deriving it. The source is a 2018 review, so the
                     reprint-route rules and the E1.1 test apply.
  B2.1 Voorhies    — B2.2 deliberately left Voorhies its own page and already proved
                     the Voorhies exponents unreachable by parallel-exponential decay,
                     so there is a real claim to test rather than restate.
Re-pick from the queue if a case's own yaml reveals a blocker. Avoid D1.1–D1.5: the
maintainer's scope decision on them is still open.

Verifier on every ready page, no exceptions. Park blocked cases with a resume: block —
never wait for me.

Republish dashboards at wind-down to the SAME URLs (never mint new ones):
  papers    https://claude.ai/code/artifact/99b52225-ce54-4487-99ec-3a420f2ac4ad  (favicon 📄)
  decisions https://claude.ai/code/artifact/237909a0-3432-4b82-8cfe-15b05019d036  (favicon ❓)

Halt = coordinator.md's checklist: integrate what's ready, update handoff,
republish dashboards, gates green, push, end.
```

---

## Why each part is there

**One case at a time** is the maintainer's instruction of 2026-08-07, and it replaced
"~5 builders live". The lean playbooks cut cost *per dispatch*, but the API limit is on
usage *rate*: N parallel builders burn N× as fast and one limit hit destroys N cases
instead of one. Serial work is not cheaper per page — it caps the blast radius, and a
case committed can never be lost.

**"Never block on me"** is standing. Park a blocked case with a `resume:` block
recording what is established, what each plausible answer changes, and which files to
touch; start the next at once. The maintainer reads the dashboards when they choose.

**Dispatch policy** is about avoiding the figure-review gate rather than servicing it
faster. Ranked: a worked example with printed intermediates beats an internal identity
beats a stated numerical result beats a digitised figure. Only the last needs a human.

**The verifier** stays on the maintainer's explicit instruction. Measured over the six
pages of 2026-08-07/08: **every single one needed fixes after verification, and on four
the verification changed the page's conclusion rather than polishing it.** It roughly
doubles cost per page and is the only reason those four are not live and wrong.

**The favicons are recorded** because a changed favicon reads as a different page — the
maintainer finds these tabs by their icon. Keep them stable.

## State at the time of writing (2026-08-08)

| | |
|---|---|
| Published pages | 64 directories, 65 catalogue entries, 73 models.yaml entries |
| Gates | `check_agreement` 64 pages / 0 regressions; `check_metadata` OK, 0 warnings |
| Queue | 65 published, 52 unclaimed, 140 needs-paper, 6 covered, 3 deferred |
| Tree | clean, everything pushed |

Counts differ between page directories and catalogue entries because `B1.1` covers
`B1.5`, `A2.1` covers `A2.2`, and `C2.10` covers `D3.4`.

Six pages were added on 2026-08-07/08: `J1.3` BET, `A3.5` Ranz–Marshall, `E1.1`
Toomey–Johnstone, `A4.1` Wilke, `A2.4` tanks-in-series, `J1.1` Langmuir. All five of
the first were rebuilt from source after a parallel wave lost them to a session limit.

## Measured cost, for planning

A complete case — builder, verifier and fixer, each a separate agent — runs about
0.9–1.1M subagent tokens and 1.5–2.5 h wall clock. Repo-wide history says ~11–12M
weighted tokens per published page all-in. Six cases in one long day. Serial pacing kept
every case inside a session window; nothing was lost after the rule changed.

Most of the repo's history (415M of 490M weighted tokens) was built on Opus 5; the
2026-08-07/08 six were Fable 5. The quality guarantee is the adversarial loop rather
than the model tier — the defects it catches are diligence-class (hidden exclusions,
false claims about a source, mislabelled independence, hard-typed numbers that drift).

## Three lessons from 2026-08-07/08 worth not relearning

1. **Rebuild from source; never inherit an abandoned build.** All five rebuilt cases
   deleted the dead wave's work, and the independent re-reads disagreed with it on
   three of five — including a digit that moved a headline.
2. **A negative claim about a source is a claim** (builder.md non-negotiable 12, added
   that day). `A3.5` said the authors "did not notice" what they print; `A2.4` said a
   chapter prints no rule it does print. Both were the page's headline contribution.
3. **Settle ambiguous glyphs by arithmetic, not pixel shape** (non-negotiable 1). And
   quote printed defects verbatim, `[sic]` and all.

## Two decisions still carried

Neither blocks work — ask once, record the answer, continue.

1. Purging three copyrighted overlays from git history (needs a force-push to a public
   repo).
2. The Elsevier API key at `~/.config/elsevier/apikey` is blocked by the Bash
   permission classifier — the highest-leverage unblock for section A.

## One item for the pymrm library, not the gallery

There is **no pure-outflow boundary condition in pymrm**. Found building `A2.4`, where
the tank count is a physical rather than a mesh parameter: zero-gradient makes
`construct_convflux_upwind` reconstruct the exit face to second order, a modelling error
worth 5.9 % at N = 2 and decaying as N^−0.87 — indistinguishable from ordinary
discretisation error unless you know to look. `A2.4` derives an exact workaround
(suppress the boundary flux with `b=1, d=0`, add outflow as a
`construct_coefficient_matrix` sink on the last cell). Worth a bc kind or a documented
recipe in pymrm itself.
