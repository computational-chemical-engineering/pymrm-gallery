# pymrm-gallery

An open-science gallery of important phenomenological models from the chemical
engineering literature, each reproduced with [`pymrm`](https://github.com/) and
shown against its original experimental data.

**Status: nine pages built, four of them validated against experiment.** The
literature survey (266 models), the Quarto scaffolding, and nine complete pages
are in place; CI executes every notebook and validates metadata.

| Page | Validation |
|---|---|
| `A4.9` Duncan–Toor ternary diffusion | **experimental** — 0.59 mole % over 28 digitised points |
| `C2.1` Xu–Froment steam reforming kinetics | **experimental** — 0.0017 in conversion over 61 digitised points |
| `B1.1`+`B1.5` Thiele modulus and Weisz–Hicks | exact/reference solutions (provenance tier 6) |
| `D2.2` Van Welsenaere–Froment runaway criteria | published reference solution — 0.054 % over 30 values (tier 6) |
| `F3.1` Hatta regimes | exact/reference solutions (provenance tier 6) |
| `F1.4` Krishna–Ellenberger large-bubble holdup | **experimental** — 13.8 % mean deviation over 63 digitised points |
| `A2.3` Taylor–Aris dispersion | 1.0e-4 vs Taylor's Eq. 25; his own capillary run to 0.04 % |
| `J1.5` LDF breakthrough | 6.6e-5 vs the exact series solution (provenance tier 6) |
| `H1.7` Wijmans–Baker solution–diffusion | **experimental** — two constants fitted to two series predict the third |

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/check_metadata.py     # metadata + provenance validation
python scripts/run_pages.py          # execute every page notebook

export QUARTO_PYTHON=$(which python) # so Quarto uses the same interpreter
quarto render                        # build the site into _site/
```

## Goal

For each catalogued model, one page that shows:

1. the model as originally published,
2. the core experimental data it was validated against,
3. a `pymrm` reproduction overlaid on that data, with a quantitative agreement
   metric,
4. what `pymrm` adds beyond the original — a relaxed assumption, a finer
   discretisation, a limit the original could not reach.

Every page is a downloadable Jupyter notebook that runs unmodified in Google
Colab.

## Who it is for

- **Researchers** who need a validated starting model for their own system and
  want to adapt rather than rebuild.
- **Students and lecturers** who want a worked, data-backed example of a named
  model.
- **Coding agents** assisting either group — see `AGENTS.md` and `models.yaml`
  (planned), which expose the catalog in machine-readable form so an agent can
  find the structurally nearest existing model and adapt it, instead of
  assembling one from API documentation.

## Planning documents

| Document | Contents |
|---|---|
| [`docs/handoff.md`](docs/handoff.md) | **Start here.** Current state, papers and API access, what to build next, and the lessons worth not relearning |
| [`docs/taxonomy.md`](docs/taxonomy.md) | Classification scheme: physical domain, mathematical structure codes (`S1`–`S13`), model tiers, build priorities, the page contract |
| [`docs/catalog-A-foundations.md`](docs/catalog-A-foundations.md) | Transport closures, catalyst particle models, reaction kinetics (98 entries) |
| [`docs/catalog-B-reactors.md`](docs/catalog-B-reactors.md) | Fixed bed, fluidised bed, bubble column, trickle bed, membrane, structured/intensified (104 entries) |
| [`docs/catalog-C-adjacent.md`](docs/catalog-C-adjacent.md) | Adsorption, crystallisation/PBE, electrochemical, biochemical, polymerisation, emerging (64 entries) |
| [`docs/pdf-requests.md`](docs/pdf-requests.md) | Papers still needed to build planned pages, in priority order, with what each one is needed for |
| [`docs/pdf-findings.md`](docs/pdf-findings.md) | What was found in the papers already supplied, and what remains to extract from each |
| [`docs/data-strategy.md`](docs/data-strategy.md) | Where experimental data comes from, the legal position on digitising figures, the provenance sidecar format |
| [`docs/blueprint.md`](docs/blueprint.md) | Repository layout, Quarto publishing, Colab compatibility, `models.yaml`, CI, licensing, contribution model |

## Scope at a glance

**266 models** catalogued across ten sections, 154 of them flagged as realistic
first-wave pages. The classification is deliberately two-dimensional:
by physical domain (what a human browses) and by mathematical structure (what an
agent matches on). The structure axis is what lets a bioengineer looking at
diffusion-limited enzyme kinetics discover that it is the same model as a
catalyst pellet.

## Immediate open questions

1. ~~**Quarto?**~~ Decided: Quarto. Rationale and its four design consequences
   are recorded in [`docs/blueprint.md`](docs/blueprint.md#0-decision-record-quarto-as-the-platform).
2. **Repository name and host** — under a personal account, a TU/e organisation,
   or a new organisation that can accept outside contributors?
3. **Scope of section J** (adjacent unit operations). Including it roughly
   doubles the audience for ~30% more work, because the structures are shared.
   It also broadens the gallery beyond "reactor modelling" — a positioning
   decision, not a technical one.
4. **First six pages** — a proposal is in
   [`docs/blueprint.md`](docs/blueprint.md#9-suggested-first-six-pages).
