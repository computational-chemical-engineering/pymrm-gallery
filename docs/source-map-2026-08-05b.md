# Source map — the 2026-08-05 third drop

Thirty-four PDFs arrived in `~/papers/pymrm-gallery/` after 12:36 on 2026-08-05,
bringing the directory to **91 files**. This is the first drop that leaves
section A: it lands in **B** (particle, deactivation, non-catalytic gas–solid),
**C** (kinetics), **D** (industrial fixed-bed cases) and touches **I**.

**Every file was opened and read from its own title page or running head.** Not
one identity in this document comes from a filename, a DOI or publisher
metadata — and in this drop that mattered five times, including once where the
publisher metadata named a *different file in the same drop*. The evidence for
each identification is in the `identified_by` field of the corresponding
`queue_cases/<ID>.yaml`.

Companion records: [`papers-on-disk.yaml`](papers-on-disk.yaml) (ID → filename,
with the warnings) and [`pdf-findings.md`](pdf-findings.md) (per-file text-layer
quality and native resolution).

Live builders on `A4.6`, `A3.6`, `A4.5`, `A2.8`, `A3.7`, `A3.8` were not touched.

---

## 1. The three candidate textbooks, answered precisely

### 1.1 `process-calculation-by-watson.pdf` — Hougen & Watson, *Chemical Process Principles*

**What it is.** Hougen, O. A. & Watson, K. M., *Chemical Process Principles*, the
**combined volume**: Part One *Material and Energy Balances* (copyright 1943),
Part Two *Thermodynamics* and Part Three ***Kinetics and Catalysis*** (both
copyright 1947). John Wiley & Sons / Chapman and Hall.

**Full book or extract? Full book, 1157 PDF pages.** Advisory board, title page,
imprint page, prefaces, table of symbols, printed Contents, all three Parts, and
the indexes. Part Three's chapters, from the printed Contents: XVIII Homogeneous
Reactions (book p. 805), XIX Catalytic Reactions (902), XX Mass and Heat Transfer
in Catalytic Beds (973), XXI Catalytic Reactor Design (1007), XXII Uncatalyzed
Heterogeneous Reactions (1049). **PDF page = book page + 16**, verified at book
pp. 913 and 941.

**Text layer usable?** *Half.* Prose is excellent and fully searchable — which is
how every negative finding below was established. **Equations are destroyed**:
eq. (69), the Thiele modulus, extracts as `s^=/(mr)=/hr\-7r`.

**Native resolution: 150 ppi, CCITT-G4 bilevel — the lowest in this
repository, half the previous floor.** Rendering at 300 or 600 is interpolation.
It is nevertheless *usable*, and this was measured rather than assumed: cropped
to the upper 45 % of book p. 941 and rendered at native 150 dpi, `−0.05952`,
`39.11`, `28.40` and `2940b = −175` are all legible, and the subscripted symbols
on p. 913 read cleanly. A **whole-page** 150 dpi render is not adequate for
subscripts. Crop first — the Chapman & Cowling lesson, one step lower.

#### What it unlocks: exactly one case, `C1.1`, and it unlocks it completely

`C1.1` (Langmuir–Hinshelwood–Hougen–Watson, **T0/P1**) is catalogued to "Hougen
& Watson (1943, 1947)". **This book is the 1947 half**, by the two people the
case is named after. Not a reprint of the result — the result.

Chapter XIX *derives* the apparatus rather than quoting it: single-site and
dual-site adsorption, the dual-site concentration built from the fraction of
unoccupied centres and the number *s* of equidistant adjacent sites (*s* = 4 for
a square lattice, 6 for a triangular one, with the factor of one half for double
counting), chemisorption with dissociation, and then the surface-reaction /
adsorption / desorption controlled rate equations in the kinetic-term ×
driving-force ÷ adsorption-group form the catalogue names. Book p. 957 works a
full example — "reaction controlled by a dual site with surface reaction in a
pure feed" — including the least-squares parameter estimation on printed data.
It footnotes its own 1943 *Ind. Eng. Chem.* 35, 529 paper, so the 1943 half of
the citation is `origin_not_consulted` and does not need requesting.

**Validation route: printed worked examples and printed data. No figure
digitisation, no maintainer.**

#### Is it one of the four the textbook-canonical class is waiting on? No.

`handoff.md` names **Bird/Stewart/Lightfoot, Taylor & Krishna, Levenspiel,
Froment & Bischoff**. Hougen & Watson is none of them. The class has **not**
moved and remains inert.

But the honest reading is more interesting than that, and it is a change from the
Chapman & Cowling verdict. Chapman & Cowling is an *origin monograph for a single
case*. Hougen & Watson is **a genuine textbook** — five chapters of Part Three
alone span C1, B1, A3, D1 and B3 — and therefore it is the *first document on
disk of the kind the class was defined around. It still unlocks one case, because
one case is all that has been checked.* Specifically, Part Three also carries:

- **Thiele's effectiveness factor** (book p. 1013), *with attribution* — "E. W.
  Thiele, Ind. Eng. Chem., 31, 916 (1939)" — the modulus defined and the η curve
  plotted. `B1.1` **is already published**, so this changes nothing there; it is
  a second independent witness if that page is ever re-audited.
- The **Gamson–Thodos–Hougen** and **Wilke–Hougen** packed-bed *j*-factor
  correlations with their *Trans. AIChE* citations — section A3 material.
- Chapter XXI *Catalytic Reactor Design* (D1 territory) and Chapter XXII
  *Uncatalyzed Heterogeneous Reactions* (B3 territory).

**Each of those is a separate `E1.1` question and none of them has been asked.**
Asking them is cheap — the prose layer is greppable — and it is the single
highest-value follow-up this drop creates. **Do not claim any of them without
running the test.** The `E1.1` failure was precisely a source that *contains* a
relation but does not name, attribute or test it, and a textbook is the document
class most likely to fail that way.

**Net: one case (`C1.1`) unlocked completely; the class stays inert; a genuine
lead worth an hour of checking.**

### 1.2 `process-calculation-by-watson (1).pdf` — a duplicate, not a paper

**Byte-identical.** Both files are 30 910 848 bytes with md5
`cad65070d6bd280902234efbe9d64f27`. It is a browser duplicate-download suffix,
exactly like the `(1)` on the `A2.7` Westerterp file. **Not** a second part, not
a second volume, not a different edition. Mapped nowhere; nothing to request.

### 1.3 `2015.205681.Industrial-Chemical.pdf` — a book, and it unlocks nothing

**What it is.** Hougen, O. A. & Watson, K. M., ***Industrial Chemical
Calculations: The Application of Physico-Chemical Principles and Data to Problems
of Industry***, **second edition**, John Wiley & Sons / Chapman & Hall;
"COPYRIGHT, 1931, 1936". A *different book by the same two authors*, and the
**predecessor** of the one above — *Chemical Process Principles*' own preface
quotes "the preface to the first edition of Industrial Chemical Calculations".

**Full book, 504 PDF pages. No text layer at all** — `pdftotext` returns nothing
and `pdffonts` lists no font. Nothing whatever could be inferred without
rasterising. **CCITT-G4 bilevel at 600 ppi native**, the highest recorded here.
Identified from its rendered title page, imprint page and printed Contents.

**Verdict: it unlocks ZERO catalogue cases,** and this is settled by its own
Contents page, not by inference:

> Part I — I Weights and Compositions; II Stoichiometry; III Ideal Behavior of
> Gases; IV Vaporization and Condensation; V Thermophysics; VI Thermochemistry at
> Standard Conditions; VII Thermochemistry of Industrial Reactions and Fuels;
> VIII Weight and Heat Balances of Combustion Processes; IX Weight and Heat
> Balances of Chemical and Metallurgical Processes.
> Part II — X Crystallization, Adsorption, and Distribution; XI Compressibility
> of Gases; XII Entropy and Free Energy; XIII Fugacity and Thermal Properties at
> High Pressures; XIV Chemical Equilibria.

It **ends at chemical equilibria** on book p. 438. **There is no kinetics chapter
and no catalysis chapter.** Nothing in the 266-case catalogue is about
stoichiometry, humidity charts or combustion heat balances, and everything it
*does* carry is superseded by Parts One and Two of the C1.1 monograph — which is
also on disk *and* has a text layer.

**Explicitly: it is not one of the four the class is waiting on either.**

### 1.4 `90131.pdf` — not a monograph at all

6 MB and a bare numeric filename, so the brief listed it as a possible book. It
is **Kiani, D. & Wachs, I. E., "The Conundrum of 'Pair Sites' in
Langmuir–Hinshelwood Reaction Kinetics in Heterogeneous Catalysis", *ACS
Catalysis* (2024), doi:10.1021/acscatal.4c02813, 11 pp, born-digital.**

It maps to no case on its own, but it is the natural modern companion to `C1.1`:
it is a critical analysis of exactly the adjacent-pair assumption Hougen & Watson
build on book p. 913. Recorded there as a secondary. **It does not unlock a
case by itself and must not be presented as the source of the rate law.**

---

## 2. What is now buildable, ranked by tier then priority

**Twenty-seven cases moved from `needs-paper` to `unclaimed`.** Ranked below by
tier first, then priority — with the validation route, because the route decides
whether a builder can finish unattended.

### Tier 0

| Rank | Case | Title | P | Validation route | Note |
|---|---|---|---|---|---|
| 1 | `C1.1` | LHHW | **P1** | **printed worked examples + printed data** | the monograph; build first — C2.3/C2.4/C2.10/C2.19 all import it |
| 2 | `C1.3` | Mars–van Krevelen | **P1** | **printed tables + the paper's own SO₂ converter calculation** | origin paper; **retires the worst api-text in the repo** |
| 3 | `C2.3` | Ammonia synthesis (Dyson–Simon) | **P1** | **printed table** (Table I constants) | pairs with `D3.2` |
| 4 | `C2.10` | o-Xylene → phthalic anhydride | **P1** | **printed constants** | pairs with `D3.4` |
| 5 | `B1.3` | Bischoff generalised modulus | **P1** | **analytical** — collapse is computable, band width measurable | completes the modulus quartet |
| 6 | `B1.7` | Mears criteria | **P1** | **derived** — sweep where each criterion actually holds | same shape as published `B1.4` |
| 7 | `B2.1` | Voorhies coking law | **P1** | printed data — **table-vs-figure unconfirmed** | head of the B2 series |
| 8 | `B2.2` | Froment–Bischoff deactivation | **P1** | **qualitative profile signature** (ascending vs descending coke) | Part II still api-text only |
| 9 | `B2.3` | Levenspiel deactivation kinetics | **P1** | **discriminability sweep** — no data needed | *not* the monograph |
| 10 | `B3.2` | Grain model | **P1** | **structural + reduces to published `B3.1`** | Part II still api-text only |
| 11 | `C1.2` | Eley–Rideal | **P1** | **analytical only — no data in the source** | see the scope caveat, §4 |
| 12 | `D3.2` | Ammonia synthesis converter | **P1** | **printed table**, bed config not in the source | scope caveat, §4 |
| 13 | `D3.4` | o-Xylene multitubular | **P1** | **printed constants**; 1D-vs-2D is the paper's own argument | best "pymrm adds something" in the drop |

### Tier 1

| Rank | Case | Title | P | Validation route | Note |
|---|---|---|---|---|---|
| 14 | `B3.4` | Sohn additive reaction times | **P1** | **analytical** — measure the approximation's own error | best-shaped page in the drop |
| 15 | `B3.9` | Broido–Shafizadeh (Di Blasi review) | **P1** | likely a **parameter compilation table** | born-digital, cheapest extraction |
| 16 | `C2.4` | Methanol synthesis (Graaf) | **P1** | **printed constants** | pairs with `D3.3`; C2.5 comparison impossible |
| 17 | `D3.3` | Methanol synthesis reactor | **P1** | **printed table + a stated "+5 %" target** | scope caveat, §4 |
| 18 | `B1.9` | Wakao–Smith random pore | P2 | **possibly EXPERIMENTAL** — 5 pellets, 1–12 atm | table-vs-figure unconfirmed |
| 19 | `B3.5` | CaO carbonation deactivation | **P1** | **printed refitted parameters** | reprint route — catalogue names a paper not on disk |
| 20 | `B3.10` | Oxygen-carrier reduction/oxidation | P2 | **printed kinetic parameters** | catalogue names the wrong journal |
| 21 | `B3.7` | Kobayashi two-competing-rate | P2 | **possibly EXPERIMENTAL** — measured yields | exponents concatenated in OCR |
| 22 | `B2.6` | Poisoning shell-progressive | P2 | **printed parameter estimates** on benzene/Ni/thiophene | reprint route — Wheeler not on disk |

### Tier 2 / 3

| Rank | Case | Title | P | Validation route | Note |
|---|---|---|---|---|---|
| 23 | `C2.19` | Ethanol dehydrogenation | **P1** | **printed tables**, two eras | catalogue had no citation at all |
| 24 | `B1.13` | Bimodal / macro–micro pellet | P2 | **derived/computational** | **proposal, not a match** — see §4 |
| 25 | `B2.4` | Beeckman–Froment pore blockage | P3 | qualitative signature | completes the B2 series |
| 26 | `B1.10` | Feng–Stewart pore network | P3 | **structural** — test the two-radius claim | downstream of published `A4.3` |
| 27 | `B3.8` | CPD / FG-DVC coal pyrolysis | P3 | **printed structure coefficients** | three papers, completely sourced |

*(Ranks 23–27 are the T2/T3 block; `C2.19` is P1 within it.)*

### The figure-route warning, again

**There is still no maintainer available for figure review as of 2026-08-05.**

Unlike the previous drop, **nothing here is figure-only outright**. But **five
cases have an unconfirmed table-vs-figure split** and a builder could reach the
review gate on any of them: `B2.1`, `B1.9`, `B3.7`, `B3.5`, `B2.4`. In every case
the check is one render away and should be done *before* dispatch, not after.

The five with an unconfirmed split are ranks **7, 18, 19, 21 and 25**. Every
other rank — 1–6, 8–17, 20, 22–24, 26–27 — needs no figure work at all.
**Dispatch those first.**

### Scope decisions to make before dispatching

1. **Four pairs share one file each.** `C2.3`+`D3.2` (Dyson & Simon),
   `C2.10`+`D3.4` (Froment 1967), `C2.4`+`D3.3` (Graaf). Build the *kinetics*
   page first and have the *reactor* page import it. **Do not dispatch two
   builders who will each transcribe the same table.**
2. **The whole B2 deactivation series is now sourced end to end** — `B2.1`
   Voorhies → `B2.2` Froment–Bischoff → `B2.3` Levenspiel → `B2.4`
   Beeckman–Froment, plus `B2.6`. That is unusual and is worth building as a
   connected ladder rather than five unrelated pages, the way `D1` is designed.
3. **`B1.9`, `B1.10` and `B1.13` are all bidisperse/network pellet transport** and
   were all unblocked together. Their scopes overlap. Decide the split.
4. **`B3.7`, `B3.8` and `B3.9` are one pyrolysis cluster** (coal two-reaction,
   coal network, biomass review) and `B3.8`'s own introduction cites `B3.7`.
   Same decision.
5. **`C1.1` before everything in C and D.** Four of this drop's cases are LHHW
   instances.

---

## 3. Catalogue-metadata errors found

Flagged the way `A1.6`, `A1.8`, `A2.9` and `A3.13` were. Fix the catalogue rows
when the pages are built; do not propagate any of these into `meta.yaml` or
`models.yaml`.

| Case | Catalogue says | Document says | Severity |
|---|---|---|---|
| `C2.5` | Vanden Bussche & Froment, ***J. Catal.* 161:1 (1996)** — with a recorded DOI `10.1002/cjce.5450740524` | that DOI is the **Canadian Journal of Chemical Engineering**. The auto-resolver returned the authors' *other* 1996 paper, and **that is the paper that was then acquired** | **wrong DOI, and it cost an acquisition** |
| `B3.10` | Abad & Adánez, ***Chem. Eng. Sci.*** (2007) | **Energy & Fuels 21(4) 1843–1853 (2007)**. Right authors, right year, wrong container | **wrong journal** |
| `D3.3` | Graaf; Vanden Bussche & Froment — DOI `10.1016/0009-2509(93)80150-o` | a **1993** CES DOI, matching *neither* document on the row (Graaf's kinetics paper is CES **43 (1988)**) | **wrong DOI**, removed |
| `B3.5` | Grasa & Abanades, IECR **45 (2006)** | on disk is Criado, Arias & Abanades, IECR **57 (2018)** 12595–12599, which *refits* the 2006 curve. The 2006 paper is absent | catalogue over-promises what one file gives |
| `B3.6` | Field **(1967)**; Baum–Street **(1971)** | neither on disk; a 1982 I. W. Smith review is, and its sufficiency is **unchecked** | stays blocked, deliberately |
| `B3.4` | *Metall. Trans.* (1978) | Metallurgical Transactions **B**, volume **9B**, March 1978, **89–96** | incomplete |
| `B2.1` | Voorhies, **A.** | the paper prints **Alexis Voorhies, Jr.** | trivial |
| `B2.3` | "**(1972)**" — bare year, no author, no journal | J. Catalysis **25 (1972) 265–272**, Levenspiel | citation too thin to be checkable; now resolved |
| `B3.7` | "**(1977)**" — bare year | 16th Symp. (Int.) Combustion (1977), Kobayashi, Howard & Sarofim, p. 411 ff. | same; now resolved |
| `C2.19` | "(in this suite's teacher exercises)" — **no literature at all** | Franckaerts & Froment, CES **19 (1964) 807–818** | case was unbuildable by AGENTS.md regardless of difficulty; now resolved |
| `C1.2` | "**—**" — no citation at all | filled from Prins, *Top. Catal.* **61 (2018) 714–721**, with a scope caveat | was unresolvable |
| `B1.13` | "**—**" — no citation at all | *proposed* Gheorghiu & Coppens, *AIChE J.* **50(4) 812–820 (2004)** | proposal, not a confirmed match |
| `C2.10` | Froment (1967); **Calderbank** | Froment 1967 on disk and correct; **Calderbank absent** | benign |
| `B2.2` | CES 16 (1961), **17 (1962)** | only the **1961 Part I** PDF arrived; Part II remains api-text only | benign |
| `B3.2` | Szekely & Evans, CES 25 (1970) | correct; **Part II (CES 26, 1971) still api-text only** | benign |

**Six DOIs were upgraded** from `crossref-auto`/`crossref-verified` to
verified-against-the-document (`B1.3`, `B1.7`, `B2.1`, `B2.4`, `B3.4`, `C2.3`,
`C2.10`), and **two were removed as wrong** (`C2.5`, `D3.3`).

### The new variant of the DOI failure, and the fix

`handoff.md` records that "a DOI resolved from a terse citation is usually
wrong", with `F1.6` and `F3.4` resolving to *book reviews*. **`C2.5` is a new
variant and a worse one**: the **year agreed**, the **authors agreed**, and the
**container did not** — and container disagreement was not being checked. The
consequence was material: the acquisition was made from the DOI, so the wrong
paper was supplied and a `T1/P1` case is still blocked.

**Add the container to the check.** Year plus authors is not enough.

---

## 4. Scope caveats that must reach the pages

- **`C1.2` (Eley–Rideal).** The source is a **historical and terminological**
  paper whose thesis is that the mechanism should be called **Langmuir**–Rideal,
  and that what Eley and Rideal actually studied is a different reaction. It
  prints the L–H rate law in full, names and distinguishes all three mechanisms,
  and attributes the third one to Langmuir. **It contains no data and no tables.**
  A page that builds the rate law and reports the attribution is defensible; one
  that presents this as the origin of the "Eley–Rideal rate law" is not. Same
  class of finding as `B1.4` ("It was shown by WEISZ [9]") and `H1.1`.
- **`D3.2`.** The document gives the rate law and the pellet correction, **not a
  multi-bed quench/intercooled converter design**. Build the converter with the
  bed configuration stated as a chosen illustrative case, and name what is
  missing — the `A1.8` / `A3.13` discipline.
- **`D3.3`.** Neither document is a **Lurgi quasi-isothermal** reactor: one is a
  spinning-basket kinetics study, the other a reversed-flow reactor. Build the
  quasi-isothermal tube from Graaf's kinetics and say no published Lurgi design
  was read.
- **`B3.5`, `B2.6`, `B3.9`** are all **reprint route** — use
  `reference_read_from` for what was read and `origin_not_consulted` for
  Grasa & Abanades (2006), Wheeler (1951) and Broido/Shafizadeh respectively.
- **`B1.13`** is a **proposal**. The paper optimises pore networks; the case
  describes two-scale resistance. The honest page is the interesting one (test
  the near-optimality of fractal hierarchies), but the choice must be made
  deliberately.
- **`B3.2`.** Szekely & Evans' own claim is that the model "reproduce[s] the
  general **trends**" using "**reasonable values**" for the parameters. That is a
  trend claim, not a fit. Do not present it as a validated prediction.

---

## 5. Things that are not what their filename, metadata or catalogue row suggests

- **`Can J Chem Eng … Vanden Bussche` is not `C2.5`, and it is not even a reprint
  route for it.** Same authors, same year, different journal, different problem —
  the `A3.5` pattern. The twist that makes it worse: the STAR paper does not use
  its own authors' kinetics. Its text says the reactor model "was combined with a
  set of kinetic equations … proposed by **Graaf** et al. (1988)", and its
  reference list carries Graaf, not *J. Catal.* 161. **It prints the competing
  model.** `C2.5` stays `needs-paper`; the file is mapped to `D3.3`.
- **`i160028a013.pdf` carried the wrong `Title` metadata** — the title of a
  *different file in the same drop*. The document is Dyson & Simon 1968. Settled
  by the **ACS download stamp printed on every page**
  (`article-pdf/7/4/605/…`), which is the cheapest reliable identity check
  available for ACS files in this repo. `ef070025k.pdf` needed it too — its
  `Title` is `No Job Name`.
- **Five files open with the previous article's text**: `B1.7` Mears, `B2.4`
  Beeckman–Froment, `C2.3` Dyson & Simon, `B3.8`'s Solomon, and the unmapped
  Kissinger. `B1.7` is the nastiest — the preceding article discusses "the
  Voorhies (1945) relationship", which makes it look like `B2.1`. **Always scroll
  past the first screen.**
- **`413a375.pdf`, a bare numeric filename, is an ordinary journal article** —
  Chem. Papers 41(3) 375–393 (1987) — with a Cyrillic abstract, not a report.
  **`90131.pdf`, also bare and numeric and 6 MB, is a 2024 ACS Catalysis paper**,
  not a monograph.
- **`BF02822675.pdf` extracts with a space between every letter**, so a keyword
  `grep` on it returns nothing and you will conclude the content is absent. It is
  not.
- **`ac60131a045.pdf` is a real, correctly-named paper that fits nothing.** See
  §6.

---

## 6. On disk, identified, mapped to no catalogue case

- **`ac60131a045.pdf`** — Kissinger, H. E., "Reaction Kinetics in Differential
  Thermal Analysis", *Anal. Chem.* **29(11) 1702–1706 (1957)**, National Bureau
  of Standards. The Kissinger method: activation energy from the shift of a
  DTA/DSC peak with heating rate. **No catalogue case corresponds to it.** It is
  *not* `C1.5` — `C1.5` is a rate-law formalism, this is a data-reduction method
  for a thermal analyser, and no case in A, B, C or J is about thermal analysis.
  **Do not force it onto `C1.5` and do not let its arrival close anything.**
  Reported, not mapped — the `ie50546a056.pdf` precedent.
- **`2015.205681.Industrial-Chemical.pdf`** — Hougen & Watson, *Industrial
  Chemical Calculations*, 2nd edn (1936). §1.3.
- **`process-calculation-by-watson (1).pdf`** — byte-identical duplicate. §1.2.
- **`90131.pdf`** — Kiani & Wachs (2024). Secondary for `C1.1`. §1.4.
- **`1-s2.0-S0009250997003850-main.pdf`** — Pan & Zhu, "Study on
  diffusion–reaction process inside a cylindrical catalyst pellet", *CES*
  **53(5) 933–946 (1998)**. A competing approximate effectiveness factor for
  arbitrary kinetics in a **cylinder**, plus Weisz–Hicks multiplicity for a
  cylinder. Secondary for `B1.3`; recorded there. `B1.2` is published and this
  does not reopen it.
- **`1-s2.0-S0082078482802816-main.pdf`** — I. W. Smith, *The Combustion Rates of
  Coal Chars: A Review* (1982). The `B3.6` near-miss; see §7.

---

## 7. What stays blocked, and why

| Case | T/P | Why still blocked |
|---|---|---|
| `C2.5` Methanol synthesis (alternative) | T1/**P1** | **needs the real paper.** *J. Catal.* 161(1) 1–10 (1996). The near-miss on disk prints the *rival* model. The recorded DOI was wrong and caused the wrong acquisition — do not repeat it from a DOI, request it **by title and container**. |
| `B3.6` Char combustion | T1/P2 | **deliberately parked, not forced.** An open question first: does I. W. Smith's 1982 review print the Field / Baum–Street correlations with attribution well enough to pass the reprint test? It carries the *subject* — intrinsic reactivities corrected for pore diffusion — and carrying the subject is exactly what `E1.1` shows is not enough. The file's text layer is the worst in the drop (37 chars on page 3) so the check needs renders. **Answer it before requesting Field (1967) and Baum & Street (1971).** This is the `A3.5` discipline. |
| `I1.7` Reverse-flow reactor | T1/**P1** | needs Matros & Bunimovich, *Catal. Rev.* **38 (1996)**. The Vanden Bussche STAR paper states the principle and simulates a reversed-flow reactor as its baseline, but is about the STAR variant on one chemistry. **Open question recorded on the case**; answer it before requesting. |
| `I1.8` Chemical looping combustion reactor | T1/P2 | needs Adánez et al., *PECS* **38 (2012)**. Two looping papers arrived and neither is a reactor model — one is particle kinetics (→ `B3.10`), one is calcium looping (→ `B3.5`). Recorded so nobody re-identifies them. |
| `B1.11` Multicomponent pellet with MS/DGM | T1/P2 | **not flipped, but its likely source arrived.** `C2.19`'s Franckaerts & Froment (1964) is ethanol dehydrogenation, which is literally what `B1.11`'s row names — but whether the paper reports pellet-scale information has **not** been checked. Check before dispatching, and do not dispatch `B1.11` and `C2.19` onto the same transcription. |

Also unchanged from the previous map: `A2.4` (needs **Levenspiel** the
*monograph* — note that Levenspiel the *author* arriving for `B2.3` does **not**
close it), `A3.2`, `A3.5`, `A4.1`, `A4.8`.

---

## 8. What to ask the maintainer for

Short, in priority order.

1. **The real `C2.5` paper — request it by title, not by DOI.**
   Vanden Bussche, K. M. & Froment, G. F., "A Steady-State Kinetic Model for
   Methanol Synthesis and the Water Gas Shift Reaction on a Commercial
   Cu/ZnO/Al₂O₃ Catalyst", ***Journal of Catalysis* 161(1) 1–10 (1996)**. The
   wrong one has already been supplied once because the recorded DOI pointed at
   the *Can. J. Chem. Eng.* paper. This is T1/P1 and it is the only blocker on
   the C2.4-vs-C2.5 comparison page the catalogue asks for.
2. **A figure-review decision, or a maintainer.** Unchanged and still the binding
   constraint on the *previous* drop (`A3.14`, half of `A3.10`). Five cases in
   *this* drop have an unconfirmed table-vs-figure split.
3. **The four textbook-canonical monographs**, unchanged: Bird/Stewart/Lightfoot,
   Taylor & Krishna, **Levenspiel** (which also closes `A2.4`), Froment &
   Bischoff. **Hougen & Watson is a real textbook and it is not one of them.**
   Note that Froment & Bischoff is now the most valuable of the four: this drop
   put **five** Froment papers on disk (`B2.2`, `B2.4`, `C2.10`/`D3.4`, `C2.19`,
   `D3.3`) and the book is `D1`'s canonical treatment.
4. **Froment & Bischoff, *CES* 17 (1962) 105** (fouling Part II) and
   **Szekely & Evans, *CES* 26 (1971) 1901** (grain model Part II). Both parts I
   arrived as scans; both parts II remain api-text only, and that route drops
   decimal points.
5. **Grasa & Abanades, *IECR* 45 (2006)** — would turn `B3.5` from a reprint
   route into an origin build.
6. Lower value: Field (1967) and Baum & Street (1971) for `B3.6` — **but only
   after the Smith-review question in §7 is answered**; Matros & Bunimovich
   (1996) for `I1.7`, same condition; Calderbank for `C2.10`; CPD part 2 for
   `B3.8`.

**Nothing in this list is urgent except item 1.** The drop moved twenty-seven
cases and the queue's binding constraint is now builder throughput and the figure
gate, not paper supply.
