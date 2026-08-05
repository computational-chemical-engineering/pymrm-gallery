# Source map — the 2026-08-05 second drop

Sixteen PDFs arrived in `~/papers/pymrm-gallery/` after 12:09 on 2026-08-05,
bringing the directory to 57 files. All sixteen are section A. This document
records what each one is, what moved, what did not, and why.

**Every file was opened and read from its own title page.** Not one identity in
this document comes from a filename, a DOI, or publisher metadata. Two of the
sixteen have no text layer at all and had to be rasterised before anything at all
was known about them. The evidence for each identification is written into the
`identified_by` field of the corresponding `queue_cases/<ID>.yaml`.

Companion records: [`papers-on-disk.yaml`](papers-on-disk.yaml) (the ID → filename
map, with warnings) and [`pdf-findings.md`](pdf-findings.md) (per-file text-layer
quality and render resolution).

---

## 1. The Chapman & Cowling question, answered precisely

**What it is.** `cambridge-mathematical-library-…-non-uniform-gases-…pdf` is
Chapman, S. & Cowling, T. G., *The Mathematical Theory of Non-Uniform Gases: An
Account of the Kinetic Theory of Viscosity, Thermal Conduction and Diffusion in
Gases*, **third edition (1970)**, prepared in co-operation with D. Burnett;
reissued in paperback in the Cambridge Mathematical Library with a Foreword by
Carlo Cercignani (1990), reprinted 1993.

**Full book or extract? Full book.** 448 PDF pages: half-title, title page,
imprint page, Cercignani's Foreword, the Prefaces to the first and third
editions, Contents, the eight-page Chapter and Section Titles list, Chapters
1–19 (book pages 1–406), the Historical Summary, name index, subject index, and
the table of references to numerical data. Nothing is missing. PDF page = book
page + 24 in the body, verified at six separate pages.

**Scan or text? Both, and the distinction matters.** ABBYY FineReader 8 over
**CCITT-G4 bilevel page images at 300 ppi native**. The prose is good and
searchable — which is how the negative findings below were established. The
numbers are not: the 1970 typesetting uses a raised mid-dot decimal separator
that the OCR renders as a dot, an apostrophe or nothing, and a leading zero is
routinely read as the letter "o". Table 22 comes back as an unreadable run of
capitals. **Render at 300 dpi, not 600** — 600 upsamples a 1-bit image.

### What it unlocks: exactly one case, `A4.6`, and it unlocks it completely

`A4.6` (Chapman–Enskog, T0/P1) had `catalog_reference: —` — no citation at all —
and had been sitting in `needs-paper` since 2026-07-31 with the note that
CrossRef could not resolve anything. This book is not a reprint of the result; it
**is** the result, by one of the two people it is named after. So `reference` is
filled and `reference_read_from` is deliberately left empty.

Asked the `AGENTS.md` question — *does it carry the result the case is for, or
does it merely mention the topic?* — section by section:

| Needed | Where | Form |
|---|---|---|
| The binary diffusivity formula | §14.2, eq. (14.2, 4); derived §9.81 | equation |
| Second approximation, Kihara form | §14.21 | equation |
| Collision integrals W₁₂⁽¹⁾(1), W₁₂⁽²⁾(2), A, B, C | **Table 6** (p. 185), **16 rows**, kT/ε₁₂ = 0.3 → **100** | **printed table** |
| Lennard-Jones ε/k, σ, rₘ | **Table 17** (12,6 model); Table 18 (exp;6) | **printed table** |
| Measured D₁₂ to validate against | **Table 22** (p. 263), ~45 gas pairs at S.T.P., each with its source reference, plus σ′₁₂ and ½(σ₁+σ₂) | **printed table** |
| D₁₂ vs temperature; self-diffusion | Tables 23, 24, 25 | **printed table** |

That is the whole build, and **every piece of it is a printed table**. Fig. 8 on
p. 185 plots the same collision integrals that Table 6 tabulates — ignore the
figure, use the table.

*A small caution, recorded because it happened during this very pass.* Table 6's
last row reads `100·0`, and read off a full-page render at 200 dpi the mid-dot
disappears and it looks like `1000`; a miscount of the rows went with it. Both
were caught only on a second, cropped read at native resolution. Everything in
this book is a mid-dot decimal — **check anything numeric twice, on a crop, at
300 dpi.** **`A4.6` therefore needs no figure digitisation and no
maintainer.** It is the single cleanest new build in this drop.

The book also supplies its own honesty caveat, which the page must carry: §14.2
opens "Observations of the mutual diffusion of pairs of gases are difficult to
make, and such observations as have been made are liable to a fairly large
experimental error", and §14.32 says the tabulated D₁₂ "may, in unfavourable
cases, be in error by several per cent."

### What it does *not* unlock — say this out loud

**It is not the textbook-canonical lever.** `handoff.md` §"Recommended next
moves" item 2 records that the ~20 T0 textbook-canonical cases are approved in
principle but inert because none of the monographs is on disk, and names the four
that matter: **Bird/Stewart/Lightfoot, Taylor & Krishna, Levenspiel, Froment &
Bischoff**. Chapman & Cowling is none of them, and neither is Li & Kwauk. That
class has not moved. What arrived is better in kind but narrower in reach: two
*origin* monographs for two specific cases, rather than a textbook whose
restatements could source many. The dashboard line "one book unlocks more than
any paper" is still an outstanding request, and the four books above are still
the ask.

**It does not close `A4.1`.** Checked deliberately, because "Wilke" appears in
the batch twice and both appearances are traps. The only Wilke in the book is
Buddenberg & Wilke (1949), cited once in §12.43 for the **viscosity** of a
mixture, with a single name-index entry at p. 244. Wilke's 1950 diffusion mixture
rule is not stated, not named and not attributed anywhere. Chapter 18 gives the
exact multicomponent (Maxwell–Stefan) coefficients and §18.51 the
one-gas-as-a-trace limit — the theory Wilke's rule approximates — but neither is
connected to Wilke. **This is the `E1.1` failure precisely**: the
equation-adjacent theory is present, the named result the case is about is not.
Building `A4.1` from Chapter 18 would mean reconstructing an attribution the
source does not make. `A4.1` stays `needs-paper`.

**Nothing else in A4.** `A4.3` (dusty gas) and `A4.4` (Knudsen/Bosanquet) contain
no pore transport of any kind and are published anyway. `A4.5` postdates the third
edition. `A4.2` is already published from Krishna & Wesselingh 1997.

**Net verdict: one case, `A4.6`, fully unlocked with a printed-table validation
route. Zero cases beyond it. The textbook-canonical class remains inert.**

---

## 2. The second monograph — Li & Kwauk, `A1.9`

`7981631866938162.pdf` — a bare sixteen-digit filename with **no text layer at
all** (`pdftotext` returns one byte, a form feed, for each of the 214 pages).
Rasterised to identify: it is Li, Jinghai & Kwauk, Mooson, *Particle-Fluid
Two-Phase Flow: The Energy-Minimization Multi-Scale Method*, Metallurgical
Industry Press, Beijing, **1994** — exactly the reference `A1.9`'s catalogue row
names. Complete book: Chapters 1–5, Acronym, Notation, References (p. 187),
Index (p. 199).

It carries the case: §2.6 "The EMMS Model" states the model by name as a
nonlinear optimisation in eight variables under nine constraints "without using
any adjustable parameters", writes out Model LG and the three slip velocities;
Chapter 3 gives the solution algorithm and stability conditions; §5.2.1 "Drag
Coefficient and Slip Velocity" is the drag content the case title asks for; and
**book p. 158 carries a printed numeric table** — (W_st)_PFC, (W_st)_FD, ε_PFC,
ε_FD, ε_ideal, ε_a, ε\* at U_g = 2.0, 3.21 and 4.0 m/s — which is a deterministic
reproduction target for the solver rather than a curve to digitise.

Two caveats on the case yaml: "EMMS drag" as a *CFD closure* (the heterogeneity-index
formulation) postdates this book, and a page must not conflate the two; and 214
pages with no searchable text is a real extraction cost.

---

## 3. What is now buildable, ranked

Twelve cases moved from `needs-paper` to `unclaimed`. Ranked by **tier first,
then priority**, with the validation route — because the route, not the tier,
decides whether a builder can finish.

### Tier 0

| Rank | Case | Title | P | Validation route | Note |
|---|---|---|---|---|---|
| 1 | `A4.5` | Fuller–Schettler–Giddings | **P1** | **printed table** — diffusion volumes; Tables II/III comparison errors | best scan in the batch; erratum still missing |
| 2 | `A4.6` | Chapman–Enskog | **P1** | **printed table** — Table 22 (~45 gas pairs), Table 6 (16 rows), Table 17 | the monograph |
| 3 | `A2.8` | Zwietering segregation | P2 | **analytical bounds** — no dataset, none needed | τ→7 OCR trap; renders only |

### Tier 1

| Rank | Case | Title | P | Validation route | Note |
|---|---|---|---|---|---|
| 4 | `A3.9` | Billet–Schultes | P2 | **printed table** — Table 2a/2b constants; Table 3 deviations as the agreement target | best-provisioned source in the drop |
| 5 | `A3.13` | Zehner–Bauer–Schlünder | P2 | **printed table** — Mou Tables 1–3 | **scope: no radiation term** (see below) |
| 6 | `A2.9` | Baldyga–Bourne engulfment | P2 | **printed tables** — both parts, Tables 1–5 each | catalogue names the wrong journal |
| 7 | `A3.11` | Dixon–Cresswell | P2 | **printed numbers** (Table 1 only) + Olbrich analytical solution | rest is figures; live builder on `A3.12` overlaps |
| 8 | `A3.10` | Rocha–Bravo–Fair | P2 | **mixed** — constants printed, holdup/Δp data **figure-only** | Part 2 (1996, mass transfer) not on disk |
| 9 | `A3.14` | Martin–Nilles | P2 | **FIGURE-ONLY** — no tables anywhere in the paper | German; will park at the review gate |

### Tier 2 / 3

| Rank | Case | Title | P | Validation route | Note |
|---|---|---|---|---|---|
| 10 | `A2.7` | Westerterp wave model | P2 | **analytical limits + structural checks** — no data in the paper, none needed | strongest "pymrm adds something" candidate here |
| 11 | `A2.10` | Compartment models from CFD | P3 | printed tables (Tables 1–2), **but see scope** | zone network comes from a CFD solve pymrm cannot do |
| 12 | `A1.9` | EMMS drag | P3 | **printed numbers** — book p. 158 worked example | the monograph; 214 pages, no text layer |

### The figure-route warning, stated once

**There is no maintainer available for figure review as of 2026-08-05.** Two of
the twelve depend on it:

- **`A3.14` (Martin–Nilles) is figure-only outright.** The string "Tabelle" does
  not occur in the ten pages; every quantitative result is in Abb. 1–11. A
  builder dispatched here will complete the extraction and stop at the review
  gate. The honest alternatives are recorded on the case yaml: build the
  two-parameter model itself and report *no* comparison against the paper's
  data, saying so plainly — or wait.
- **`A3.10` (Rocha–Bravo–Fair) is half figure-only.** Tables I–IV carry the
  packing constants and the experimental conditions; the holdup and pressure-drop
  measurements themselves are in Figures 3–13. A page that exercises the
  correlation on the printed conditions is fine; a page that validates against
  the measurements is not, without a maintainer.

Ranks 1–7 and 10–12 need no figure work at all. **Dispatch those first.**

### Two scope decisions to make before dispatching

1. **`A4.5` + `A4.6` are one page, not two.** Chapman & Cowling's Table 22 gives
   measured D₁₂ for ~45 gas pairs; Fuller et al. give a competing correlation and
   a printed comparison against eight others, including the
   Hirschfelder–Bird–Spotz engineering form of the very Chapman–Enskog result
   Table 6 tabulates. One measured-D₁₂ axis with both predictions overlaid is the
   `A1` "one page, four correlations, one dataset" pattern the catalogue itself
   recommends, and it is the strongest single opportunity this drop creates.
   Decide before dispatching two builders who will each digitise half of it.
2. **`A3.9` + `A3.10` overlap on arranged packings.** Billet–Schultes covers both
   dumped and arranged with full constant tables and aggregate deviations but no
   raw data; Rocha–Bravo–Fair covers structured packings with constants printed
   and data in figures. Their strengths are complementary and their scope
   overlaps. Same decision.

### The `A3.13` scope caveat — this is `A1.8` repeating

`A3.13`'s catalogue line reads "Zehner–Bauer–Schlünder … Stagnant bed
conductivity **incl. radiation**". **Neither document on disk delivers
radiation.**

- Kandula (2010) is about Zehner & **Schlünder** (1970), the point-contact ZS
  model. Radiation and contact conduction are Bauer's 1978 additions.
- Mou et al. (2025) print, as their Eq. (7), the **simplified** ZBS (S-ZBS) model
  with the Smoluchowski correction folded in. Their own Introduction states that
  Syamlal & Gidaspow simplified ZBS "by ignoring the Smoluchowski **and
  radiation** effects" to obtain S-ZBS. Setting k_g\* → k_g in Eq. (7) recovers
  S-ZBS, so the model is complete and evaluable — just not radiative.
- The full ZBS with radiation is in Mou's **Supplementary Materials** (Eq. S2,
  Table S2). **The supplement is not on disk** — the PDF is the 17-page main
  article. The article refers to supplementary items nine times and none of them
  is in the file. Do not cite Eq. S2 or Tables S1–S9 as read.

Build the ZS/S-ZBS core against Mou's printed tables, name the absent component
on the page, and name the document that would complete it — exactly as `A1.8` did
with the missing Gidaspow and Wen–Yu closures. The completing document is Bauer &
Schlünder (1978), *Int. Chem. Eng.* 18, 189, or the German original Zehner &
Schlünder, *Chem. Ing. Tech.* 42 (1970) 933–941.

---

## 4. Catalogue-metadata errors found

Flagged the way `A1.6` and `A1.8` were flagged. Fix the catalogue rows when the
pages are built; do not propagate any of these into `meta.yaml` or `models.yaml`.

| Case | Catalogue says | Document says | Severity |
|---|---|---|---|
| `A2.9` | Baldyga–Bourne, ***Chem. Eng. Sci.*** (1989) | **The Chemical Engineering Journal 42 (1989)** 83–92 and 93–101, ISSN 0300-9467 — a different Elsevier journal | **wrong journal** |
| `A2.8` | *Chem. Eng. Sci.* 11 (1959) | correct, but the article runs to a printed folio **15**, not the universally cited 11. `11(1) 1–15` | incomplete |
| `A3.10` | IECR 32 (1993), **35 (1996)** | only the 1993 Part 1 (hydraulics) is on disk; **Part 2 is the mass-transfer half** and is absent — a page built from this file does not cover "hydraulics/MT" | catalogue over-promises what one file gives |
| `A3.9` | *Chem. Eng. Res. Des.* (**1993**, 1999) | Trans IChemE 77A (1999) 498–504 is on disk and is correct (Trans IChemE Part A *is* CERD); the 1993 half is absent but is superseded by the 1999 "Updated Summary" | benign |
| `A3.13` | "(1970s)" | two distinct origin papers — Zehner & Schlünder 1970 (Chem. Ing. Tech. 42, 933–941, German) and Bauer & Schlünder 1978 (Int. Chem. Eng. 18, 189). Only the latter has the radiation the case title claims | citation too vague to be checkable |
| `A2.10` | Bezzo & Macchietto (2004) | correct — two authors, Part II. **Part I is a three-author paper** (Bezzo, Macchietto & Pantelides) and is not on disk | benign, but say which was read |
| `A4.6` | "—" | no citation at all was recorded; now filled from the book | was unresolvable, now resolved |

**Three DOIs were upgraded from `crossref-auto` to verified-against-the-document**
(`A2.7` 10.1002/aic.690410902, `A3.11` 10.1002/aic.690250413, `A2.10`
10.1016/j.compchemeng.2003.08.010) and one from `crossref-verified` to
verified-against-the-document (`A4.5` 10.1021/ie50677a007). All four were right.
That is worth recording alongside the `F1.6`/`F3.4` book-review failures: the
auto-resolver is not unreliable, it is unreliable *when `catalog_reference` names
a monograph*.

---

## 5. Things that are not what their filename suggests

The four Priority-3 checks, applied to every file.

- **`Literature-PredictionofMassTransferColumns.pdf` is the paper**, not a
  literature review of it. The `Literature-` prefix is the same shape as the
  `A3.5` block (a one-page Citation Classic reminiscence and a companion study on
  a different problem), so it was rendered and read: the by-line is Billet and
  Schultes themselves, and the subtitle "Updated Summary of the Calculation
  Method of Billet and Schultes" is the authors describing their own method.
- **`20100036467.pdf` is a real report.** A NASA-style accession number was
  checked deliberately. It is a genuine 15-page technical report (its own Title
  metadata is `KSC-2010-007.pdf`) with an abstract, nomenclature, five figures
  and a reference list — not an NTRS abstract record and not an entitlement
  preview of the kind that fooled an agent on `B1.2`.
- **`7981631866938162.pdf`, a bare numeric filename, is a 214-page book.** With
  no text layer, nothing whatever could be inferred without rasterising.
- **`ie50546a056.pdf` is a Wilke paper, and it is the WRONG Wilke paper.** It is
  Wilke & **Lee**, *Ind. Eng. Chem.* 47(6) 1253–1257 (1955) — five years after,
  different journal, different co-author, different equation. `A4.1` is Wilke,
  *Chem. Eng. Prog.* 46 (1950), the multicomponent mixture rule. **This file maps
  to no catalogue case.** Its value is as a secondary source for `A4.5`/`A4.6`
  (Fuller et al. benchmark against it by name). It has not been allowed to close
  `A4.1`.
- **`1-s2.0-S0009250912007099-main.pdf` is a Stefan-tube paper that does not
  serve `A4.8`.** Mills & Chang (2013) is a *binary* 2-D Navier–Stokes +
  Maxwell–Stefan computation written to rebut the Kerkhof–Geboers equations. It
  contains no experimental data and no tables of any kind. `A4.8` is catalogued
  as the *classic MS validation experiment*, which needs measurements — the usual
  source is Carty & Schrodt (1975). `A4.8` stays `needs-paper`, with the 2013
  paper recorded so nobody reads it twice.
- **Two files open with the previous article's text.** Dixon & Cresswell's PDF
  page 1 begins with the reference list of an LNG-density paper; Wilke & Lee's
  begins with the tail of a soil-stabilisation paper. Both are complete and
  correct; a reader who checks only the first screen will reject them.
- **The `(1)` in the Westerterp filename is a browser duplicate-download suffix,
  not a part number.** There is one copy of that paper on disk.
- **The Mou filename decodes to 2024, the article is 2025.** PII
  `S0032591024012488` carries the article-in-press year; the issue is *Powder
  Technology* **453 (2025)** 120604. **Publication status checked** per the
  published-work-only rule: issue-assigned with a DOI, not a preprint, not in
  press. Kandula is a released NASA agency report — the same category as the
  DOE/METC report `A1.8` is built from and that was accepted there as the
  citable published source.

---

## 6. What remains blocked in section A, and why

Counted from `queue_cases/A*.yaml` after this pass, section A stands at **14
published, 4 covered, 4 `ready` (live builders), 12 unclaimed, 9 needs-paper** —
43 cases. Before this pass it was 21 needs-paper; **twelve moved**.

| Case | T/P | Why still blocked |
|---|---|---|
| `A1.2` Kozeny–Carman | T0/P1 | `covered` — not blocked |
| `A1.3` Darcy–Forchheimer | T0/P1 | `covered` |
| `A1.4` Eisfeld–Schnitzlein | T2/P2 | `covered` |
| `A2.2` Wehner–Wilhelm | T0/P1 | `covered` |
| `A2.4` Tanks-in-series | T0/P1 | needs-paper; `catalog_reference` is "Levenspiel" — a **monograph**, and one of the four the textbook-canonical class is waiting on |
| `A3.2` Higbie penetration | T0/P1 | needs-paper; *Trans. AIChE* 31 (1935), pre-DOI, no route found |
| `A3.5` Ranz–Marshall / Frössling | T0/P1 | **not unblocked and still not checked.** The 2026-08-05 first drop left an open question nobody has answered: whether Charlesworth & Marshall (1960, on disk) reprints enough of the correlation to source it by the reprint route. Ask that before requesting *Chem. Eng. Progr.* 48, 141–146 / 173–180 |
| `A4.1` Fick + Wilke mixture rule | T0/P1 | needs-paper. **Two near-misses arrived in this drop and neither serves it** — see §1 and §5. Remaining reprint candidate on disk: Krishna & Wesselingh 1997 (`A4.2`'s source); read its treatment of the effective-diffusivity approximation before requesting Wilke 1950 |
| `A4.8` Stefan tube / Arnold | T0/P1 | needs-paper. Mills & Chang 2013 arrived and does not serve it; needs a paper with measured ternary compositions |
| `A3.6`, `A3.7`, `A3.8` | T0–T1/P1 | **stale `needs-paper` — not a real blocker.** Their papers arrived in the 2026-08-05 *first* drop and are recorded in `papers-on-disk.yaml` (Calderbank & Moo-Young; van 't Riet; Onda et al.), but the three case files were never flipped out of `needs-paper`. No builder is live on them. **Someone should verify and flip these three** — this pass did not, to avoid colliding with whoever is mid-dispatch |
| `A1.1`, `A1.5`–`A1.8`, `A2.1`, `A2.3`, `A2.5`, `A3.4`, `A4.2`, `A4.3`, `A4.4`, `A4.7`, `A4.9` | — | **published** (14) |
| `A2.6`, `A3.1`, `A3.3`, `A3.12` (`ready`) and `A3.15` (`needs-paper`) | — | **live builders as of 2026-08-05 — not touched by this pass** |

The residue is small and its shape is now clear. Of the nine `needs-paper` cases,
one (`A3.15`) has a live builder and three (`A3.6`, `A3.7`, `A3.8`) are stale
labels on cases whose papers are already on disk. **Five section-A cases are
genuinely paper-blocked** (`A2.4`, `A3.2`, `A3.5`, `A4.1`, `A4.8`), and of those,
one (`A2.4`) needs a *monograph* that is already on the outstanding request list,
one (`A3.5`) needs a question answered before anything is requested at all, and
one (`A4.1`) has an unexhausted reprint route on disk. **Only `A3.2` and `A4.8`
are clean "please supply a PDF" asks.**

## 7. What to ask the maintainer for

In priority order, and short — the previous request lists were long because
section A was mostly blocked. It is not any more.

1. **A figure-review decision, or a maintainer.** Two unblocked cases (`A3.14`,
   and half of `A3.10`) will stall at the review gate. This is now the binding
   constraint on a P2/T1 case, not a paper shortage.
2. **The four textbook-canonical monographs**, unchanged from `handoff.md`:
   Bird/Stewart/Lightfoot, Taylor & Krishna, **Levenspiel** (which would also
   close `A2.4`), Froment & Bischoff. The two books that arrived are origin
   monographs for two specific cases, not textbooks — the class is still inert.
3. **Fuller et al.'s erratum**, `10.1021/ie50680a601`, *Ind. Eng. Chem.* 58
   (1966). It is one page, it is at the same publisher as a paper just supplied,
   and it corrects exactly the atomic diffusion volumes that are the whole
   content of `A4.5`.
4. **Carty & Schrodt (1975)**, *Ind. Eng. Chem. Fundam.* 14, 276–278 — the only
   thing standing between `A4.8` (T0/P1) and a build.
5. **Higbie (1935)**, *Trans. AIChE* 31 — `A3.2`, T0/P1.
6. Lower value, only if convenient: Rocha–Bravo–Fair Part 2 (*IECR* 35 (1996)
   1660–1667) to complete `A3.10`; Bauer & Schlünder (1978) to give `A3.13` its
   radiation term; Bezzo, Macchietto & Pantelides Part I for `A2.10`.
