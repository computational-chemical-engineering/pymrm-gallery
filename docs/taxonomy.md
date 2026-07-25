# Taxonomy and structuring scheme

This document defines *how* the gallery is organised, before listing *what* goes
in it. Three orthogonal classification axes are used, plus a priority scheme.

---

## Axis 1 — Physical domain (the browsing axis)

This is what a visitor clicks on. Ten top-level sections, each with subsections:

| # | Section | Scope |
|---|---|---|
| A | Transport & closure relations | Dispersion, mass/heat transfer coefficients, pressure drop, effective properties, diffusion models |
| B | Catalyst particle & intraparticle transport | Thiele/effectiveness, pore models, deactivation, non-catalytic conversion |
| C | Reaction kinetics | LHHW, microkinetics, named industrial rate laws |
| D | Fixed / packed bed reactors | 1D–2D, pseudo-homogeneous → heterogeneous, runaway, dynamics |
| E | Fluidised beds | Bubbling, turbulent, circulating, spouted |
| F | Gas–liquid & slurry reactors | Bubble columns, stirred tanks, airlift, slurry FT |
| G | Three-phase packed beds | Trickle bed, packed bubble, monolith slurry |
| H | Membrane reactors & separations | Pd, zeolite, polymeric, perovskite, membrane-assisted |
| I | Structured & intensified reactors | Monolith, foam, micro, reverse-flow, chemical looping |
| J | Adjacent unit operations | Adsorption/PSA, absorption/reactive separation, crystallisation/PBE, electrochemical, bio, polymer |

## Axis 2 — Mathematical structure (the pymrm axis)

**This is the axis that makes the gallery useful to a coding agent.** Every entry
is tagged with one or more structure codes. An agent asked to "build me a model
of X" first identifies the structure, then copies the nearest gallery page.

| Code | Structure | pymrm ingredients |
|---|---|---|
| `S1` | ODE-IVP in time (batch, CSTR dynamics) | `solve_ivp`, or backward Euler + `newton` |
| `S2` | ODE-IVP in space (ideal PFR marching) | marching loop + `newton` per step |
| `S3` | 1D steady BVP (particle, film) | `construct_grad` + `construct_div` + `spsolve`/`newton` |
| `S4` | 1D transient PDE (axial dispersion, breakthrough) | + accumulation matrix, time loop |
| `S5` | 1D convection-dominated (sharp fronts) | `construct_convflux_upwind` + `interp_cntr_to_stagg_tvd` deferred correction |
| `S6` | 2D PDE (radial dispersion, Graetz, monolith) | per-axis operators on `(n_z, n_r, n_c)`, `nu=1` radial |
| `S7` | Multi-domain / multi-phase coupling | `(n_z, n_phase, n_c)` layout, or `update_csc_array_indices` monolithic assembly |
| `S8` | Nested-scale coupling (reactor ↔ particle) | Schur complement via `splu`, `shapes_d` boundary unknowns |
| `S9` | Implicit multicomponent flux (Maxwell–Stefan, dusty gas) | per-face linear solve for fluxes, `NumJac(axes_blocks=[-1])` |
| `S10` | DAE / constrained (pressure–velocity, electroneutrality, equilibrium) | monolithic Jacobian with algebraic rows |
| `S11` | Population balance (internal coordinate) | treat size/age as an extra axis; growth = convection, so `S5` machinery applies |
| `S12` | Moving boundary / shrinking core | front tracking, or fixed-grid with a phase indicator |
| `S13` | Non-standard geometry | `construct_div(nu=callable)` for arbitrary area profiles |

Structure codes cluster the catalog in a way the physical taxonomy does not: a
chromatography breakthrough curve (`S4`+`S5`) and a fixed-bed tracer RTD
(`S4`+`S5`) are the *same* pymrm model with different closures. Cross-linking on
`S` codes is what lets a researcher in one field reuse a page from another.

## Axis 3 — Model rank (what "important" means here)

Rather than fake precise citation counts, entries are ranked in four tiers by
their role in the literature. Tier is a judgement call and should be revisited
per entry when the page is written.

- **T0 — Canonical.** Named after its authors, taught in every reaction
  engineering course, still cited as the reference formulation decades later.
  (Thiele modulus, Ergun, Danckwerts BCs, Kunii–Levenspiel, Monod.)
- **T1 — Field standard.** The default model in its subfield; a paper in that
  subfield is expected to use it or justify not using it. (Xu–Froment,
  Doyle–Fuller–Newman, ASM1, Wakao–Funazkri.)
- **T2 — Established alternative.** Well-cited, competes with the T1 model,
  worth showing as a comparison. (Numaguchi–Kikuchi vs Xu–Froment;
  Syamlal–O'Brien vs Gidaspow drag.)
- **T3 — Specialised / emerging.** High quality, narrower use, or recent.
  (Electrified reforming, CO2 electroreduction models, EMMS.)

## Priority scheme for building pages

Each entry carries a build priority independent of tier, because "important" and
"good first gallery page" are different questions:

- **P1** — build first. High tier, clean published data, maps to `S1`–`S6`,
  reproducible in < 200 lines. These establish the gallery's credibility.
- **P2** — build second. Needs `S7`–`S10`, or data must be digitised from figures.
- **P3** — aspirational. Needs data that may not exist in usable form, or a
  structure pymrm does not yet cover well (`S11`, `S12`, full CFD).

## Page contract

Every gallery page, regardless of section, shows the same four things in the
same order. This uniformity is what makes the gallery scrapeable by an agent.

1. **The published model** — governing equations as stated in the original,
   with the original symbols, and a table mapping them to code variables.
2. **The data** — the digitised experimental points or the original tabulated
   values, in a plain CSV committed next to the notebook.
3. **The pymrm reproduction** — one plot overlaying pymrm output on the data,
   plus a quantitative agreement metric (RMSE / parity plot / % deviation).
4. **What pymrm adds** — the improvement over the original: relaxed assumption,
   finer discretisation, a limit the original could not reach, or a sensitivity
   the original did not report.

Item 4 is the reason the gallery is more than a reimplementation archive, and it
is also where honesty matters most: if pymrm merely *reproduces* and adds
nothing, the page should say so.

See [`data-strategy.md`](data-strategy.md) for how item 2 is sourced — it is the
binding constraint on the whole project.
