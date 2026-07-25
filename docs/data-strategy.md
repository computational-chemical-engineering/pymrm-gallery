# Data strategy — the binding constraint

The gallery's promise is "core experimental data + pymrm result on one page."
Sourcing that data legally and reproducibly is harder than writing the models,
and it determines which pages get built first. This document is the policy.

---

## 1. Source hierarchy

Prefer sources in this order. Higher tiers cost less effort *and* less legal
ambiguity.

**Tier 1 — Open repositories with explicit licences.** No ambiguity, no
digitising, redistributable.
- 4TU.ResearchData (the group's own `particle_model` deposit already lives here)
- Zenodo, Figshare, Dryad
- **NETL/PSRI CFD Challenge Problems** — verified: experimental data for a
  0.305 m × 16 m CFB riser and a 0.92 m × 9 m bubbling bed, released publicly at
  `mfix.netl.doe.gov/challenge`. The best open dataset in the catalog for
  section E.
- NIST (thermophysical properties, kinetics database)
- Data-journal papers (*Data in Brief*, *Scientific Data*) — e.g. the open
  shallow-packed-bed flow/heat-transfer dataset found during the survey
  (9 experiments, CSV, ~60 Hz logging) fits section A3.

**Tier 2 — Tables printed in papers.** Numbers already in numeric form. Legally
the cleanest non-repository route (they are plainly facts) and no digitising
error. Kinetics papers (catalog C2) are rich in these — a rate table or a parity
plot's underlying values.

**Tier 3 — Supplementary material.** Increasingly common, often CSV/XLSX. Check
its licence separately from the article; SI frequently inherits the article
licence.

**Tier 4 — Digitised from figures.** The fallback for classic pre-2000 papers,
which is most of the T0 catalog. See §2.

**Tier 5 — Author contact.** For an open-science gallery that cites and credits
them, many authors will supply original data and some will co-author the page.
Worth doing for flagship pages. Slow, so start early.

**Tier 6 — Synthetic/reference-solution "data".** When no experimental data
exists, validate against an analytical solution or a high-fidelity reference
computation instead. **This must be labelled as such on the page** — never
presented as experimental. Several A-section entries will be in this category.

## 2. Digitising figures — legal position

Not legal advice; confirm with the TU/e library / research data officer, who
handles exactly this question routinely. That said, the position is well
established and favourable:

**The data points are facts, and facts are not copyrightable.** In the US this
is settled by *Feist Publications v. Rural Telephone* (1991), which rejected the
"sweat of the brow" doctrine — scientific measurements are uncopyrightable
regardless of the effort spent collecting them. The *figure* (the image, its
styling, layout, and any creative arrangement) is copyrighted; the underlying
coordinates are not.

**So: re-plot, never re-publish the image.** Extract the points, commit them as
CSV, and draw your own figure. Do not paste the original figure image into the
gallery, and do not redistribute the PDF. This single rule resolves most of the
risk.

**EU adds a database right, and an exception that covers you.** The EU *sui
generis* database right (Directive 96/9/EC) protects substantial investment in a
database independently of copyright, so the US "facts are free" analysis is not
by itself sufficient in the Netherlands. However, the DSM Directive (EU)
2019/790 Art. 3 creates a mandatory text-and-data-mining exception for research
organisations, covering both reproduction and extraction from databases, for
scientific research, on lawfully accessible works — and it **cannot be overridden
by contract**, so publisher terms of use do not remove it.

The Netherlands implemented this as **Art. 15n Auteurswet**, in force since
7 June 2021. It applies to research organisations, requires lawful access
(institutional subscription counts), permits retention of the copies with
appropriate security, and is not subject to compensation. TU/e is squarely
within scope. Art. 15o is the narrower general-purpose TDM exception with an
opt-out; you should not need it.

**Practical rules for the gallery:**

1. Only digitise from papers you have lawful access to. Institutional
   subscription or open access — never a pirated copy. This is the one condition
   the exception actually turns on.
2. Re-plot from extracted coordinates; never reproduce the source figure image.
3. Cite precisely: "digitised from Fig. 4 of Xu & Froment (1989), *AIChE J*
   35:88, DOI …". Attribution is both an academic obligation and, under the
   EU research exception, good practice.
4. Extract what the page needs, not the paper's entire data corpus. A handful of
   curves is comfortably clear of "substantial part of a database"; a systematic
   harvest of everything a publisher has ever printed is a different question.
5. Record provenance in a machine-readable sidecar (§3) so any later licence
   question can be answered per-dataset without re-tracing the work.
6. Prefer Tier 1–3 whenever a choice exists, so digitised data is a minority of
   the gallery rather than its backbone.
7. For genuinely borderline cases — a dataset that *is* the commercial product,
   e.g. DIPPR or a proprietary correlation database — do not include it; link to
   it and reproduce the model against something else.

Open-access articles are simpler still, but check the flavour: **CC-BY** allows
redistribution with attribution outright. **CC-BY-NC** or **CC-BY-ND** restrict
derivatives and commercial use — but since you are extracting *facts* and
re-plotting rather than redistributing their expression, the analysis above
still applies. Record the licence either way.

**The one thing to avoid:** bulk-mirroring PDFs into the repository. Keep the
repo free of source PDFs; store only DOIs, extracted CSVs, and provenance.

## 3. Provenance sidecar (required for every dataset)

Each `data/*.csv` gets a `data/*.meta.yaml` next to it:

```yaml
dataset_id: xu-froment-1989-fig4
title: SMR rate vs conversion, Ni/MgAl2O4, 848 K
source:
  authors: Xu, J.; Froment, G. F.
  year: 1989
  container: AIChE Journal 35(1) 88-96
  doi: 10.1002/aic.690350109
  access: institutional subscription (lawful access)
  licence: all rights reserved (facts extracted; figure not reproduced)
acquisition:
  method: digitised            # digitised | table | supplementary | repository | author | synthetic
  tool: WebPlotDigitizer 4.7
  figure: Figure 4
  operator: <name>
  date: 2026-08-01
  estimated_error: "±2% on rate axis (log scale, gridline calibration)"
columns:
  - {name: conversion, unit: "-",           description: CH4 conversion}
  - {name: rate,       unit: mol/(kg_cat s), description: observed rate}
redistribution_basis: >
  Numerical values are facts, not protected by copyright (Feist, 1991);
  extraction for scientific research is covered by Art. 15n Auteurswet
  (DSM Art. 3). Original figure image is NOT reproduced.
```

This is also the file an agent reads to decide whether it may reuse the dataset,
which makes the whole gallery machine-auditable.

## 4. Honest assessment of feasibility per section

| Section | Data outlook | Notes |
|---|---|---|
| A — transport closures | **Good** | Correlations published with the fitting data; often tabulated |
| B — particle | **Mixed** | η vs φ curves are usually theoretical, not measured; expect Tier 6 for many |
| C — kinetics | **Very good** | Rate tables and parity plots are the norm |
| D — fixed bed | **Good** | Classic hot-spot cases have printed axial profiles |
| E — fluidised bed | **Good** | NETL/PSRI challenge data is open and substantial |
| F — bubble column | **Very good** | Holdup/k_L a datasets are abundant and often tabulated |
| G — trickle bed | **Mixed** | Δp/holdup data plentiful; reaction data sparser |
| H — membrane | **Good** | Permeation data commonly tabulated; in-house data available |
| I — structured | **Mixed** | Converter light-off curves digitisable; foam data sparse |
| J — adjacent | **Good to excellent** | Breakthrough curves, battery discharge curves, ASM/ADM benchmark simulations are all well-published; battery data especially is openly available |

**Conclusion.** Data availability is good enough that the gallery is viable, but
it should be *sequenced by data availability rather than by model importance*.
The first wave should be C2 (kinetics), F1 (bubble column correlations), A4.9
(Duncan–Toor), and the in-house H1.12 / B1.12 pages — all Tier 1–3 data. Pages
depending on heavy digitisation come later, once the extraction workflow and the
provenance sidecar have been proven on a few examples.
