# Case queue

One YAML per catalogued case, tracking it from `unclaimed` to `published`.

This exists so several agents can work cases **in parallel**. Each agent touches
only its own `<ID>.yaml` and its own `<ID>/page/` staging directory, never
`models.yaml` or anything under `pages/`, so there is no race on shared files or
the git index. An integrator merges finished pages in.

```
python scripts/case_queue.py seed      # create entries for catalogued IDs
python scripts/case_queue.py status    # counts by status
python scripts/dashboards.py           # regenerate both maintainer dashboards
```

## Statuses

| status | meaning |
|---|---|
| `unclaimed` | nobody has looked at it yet |
| `in-progress` | an agent is working it now |
| `needs-paper` | the source is neither on disk nor open access — surfaces on the papers dashboard |
| `needs-input` | a decision or visual check is owed by the maintainer — surfaces on the input dashboard |
| `ready` | built and green, waiting to be merged into `pages/` |
| `published` | live on the site |
| `deferred` | deliberately not built; `blocker.detail` says why |

## Why agents halt rather than guess

A case that needs a judgement call is **halted, not abandoned**: the agent writes
its question (and, for a figure extraction, a PNG overlay under `<ID>/review/`)
into the queue entry and stops, and the next case starts. Nothing is stuck
silently, and nothing is guessed to keep a pipeline moving — which matters here
more than throughput, because a fabricated number is worse than a missing page.
