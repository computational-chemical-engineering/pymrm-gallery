# Source sweep, 2026-08-02

The queue had 266 cases: 36 published, 215 `needs-paper`, 5 `unclaimed`. All five
of those five carried `doi_source: crossref-auto` — auto-resolved, never checked
against a document. This sweep went after routes `find_papers.py` had never
tried: government repositories, institutional repositories, and the authorised
Elsevier PII full-text API.

**Headline.** Three cases are now genuinely buildable with a verified source in
hand. Two of the five "unclaimed" were never buildable and have been moved to
where they belong. And **two more crossref-auto DOIs turned out to be wrong** —
one pointing at a paper with none of the case's physics in it, two pointing at
one-page book reviews of the monographs they were supposed to cite.

The one thing that did *not* work is the thing that looked most promising: the
Elsevier API reaches the text of dozens of blocked classics but **not their page
images**, and on every scanned original the text loses its decimal points. No
case is unblocked by it. That is written up in
[`pdf-findings.md`](pdf-findings.md#the-elsevier-api-text-is-not-safe-after-1980--measured-on-15-articles-2026-08-02).

---

## 1. Buildable now — 3 cases

| Case | Source | How identity was verified |
|---|---|---|
| `A1.8` | Syamlal, Rogers & O'Brien, *MFIX Documentation Theory Guide*, DOE/METC-94/1004, Dec 1993, 54 pp. `10.2172/10145548` | PDF downloaded from OSTI and opened; its own title page names all three authors, the report number and the date. |
| `A4.7` | Krishna & Baur, *Modelling issues in zeolite based separation processes*, Sep. Purif. Technol. **33**(3) 213–254 (2003). `10.1016/S1383-5866(03)00008-X` | Elsevier full text; header gives title, journal, volume, issue, pages, date; body gives the author list. |
| `H1.9` | same paper as `A4.7` | same |

### `A1.8` — a real PDF, but it carries one drag law of three

OSTI is a US DOE public repository and the MFIX Theory Guide is a government
technical report, so this was a clean download. Section 2.2.1 gives the
Syamlal–O'Brien law complete: the terminal-velocity-to-drag conversion credited
to Syamlal & O'Brien (1987) as eq. (11), Garside & Al-Dibouni's closed-form
*V*ᵣₘ as eq. (12) with *A*, *B* in (13)–(14), and Dalla Valle's single-sphere
*C*_D.

It does **not** contain Wen–Yu — the strings "Wen and Yu" and "Wen & Yu" do not
appear anywhere in it — and cites Gidaspow only in passing (Gidaspow 1986, for
needing a low-solids-fraction correlation). So the case's title over-promises
against this one document. In practice the gap is small: Wen & Yu (1966) is
already on disk as the `A1.6` file and Ergun (1952) as the `A1.1` file, so two of
the three are covered from papers the gallery holds. Only Gidaspow's *blend rule*
(Ergun below ε ≈ 0.8, Wen–Yu above) still needs the 1994 monograph — and it must
not be written from memory.

### `A4.7` and `H1.9` — one PII, and it sourced neither

The brief's question was whether one 1995 paper really sources both cases. It
sources **neither**, and the mapping was wrong in a way that matters.

`10.1016/0009-2509(95)00102-B` is van den Broeke & Krishna, *Experimental
verification of the Maxwell–Stefan theory for micropore diffusion*, CES **50**(16)
2507–2522. (Note the author order: the catalogue writes "Krishna & van den
Broeke"; van den Broeke is first author.) The paper reports **breakthrough
experiments in a packed bed** of Kureha microporous activated carbon and a Takeda
carbon molecular sieve.

- For `H1.9` (*Maxwell–Stefan zeolite membrane*) it is simply the wrong paper.
  There is no membrane in it and no permeation model. Every occurrence of
  "membrane", "zeolite" and "silicalite" in its full text sits in the
  introduction or the reference list — Kapteijn/Bakker/Zheng on silicalite
  membranes, Rao & Sircar on carbon membranes. This is the `F3.5` / `G1.7`
  failure mode again: a DOI that resolves cleanly to a document that does not
  carry the case.
- For `A4.7` (*Krishna zeolite/micropore MS*) it is closer — the physics is
  micropore Maxwell–Stefan diffusion — but the adsorbents are carbons, not
  zeolites, and its text has a specific hole: **Table 4, the single-component
  Maxwell–Stefan diffusivities, is absent from the API text altogether**, present
  only as five cross-references from the body. Tables 1–3 do come through with
  decimal points intact.

Both cases are now pointed at the catalogue's *other* named reference, Krishna &
Baur (2003), which does carry them: its stated second objective is "to show how
the M–S approach can be incorporated into the modelling and design of practical
devices such as membrane permeation units and fixed bed adsorbers", and it works
MFI-membrane separation of hexane isomers explicitly.

**On text quality, since the brief asked specifically.** 1995 being
"post-scan-era" is not true for Elsevier — CES was still being scanned, and the
1995 text shows it: species labels corrupted (`CaH6` for C₃H₆, `C3Hs` for C₃H₈,
`COz` for CO₂), the corrected diffusivity Ð̄ rendered `/)~(0)`, and a whole table
missing. Decimal points *do* survive in it, which is better than the 1970 papers,
but that is not enough to call it clean. Krishna & Baur (2003) is a different
animal: born-digital, with `Ð`, `Θ`, `θ`, `Γ`, `Pa−1` and five tables arriving
with their headers intact. That one is usable for numbers, and it is why both
cases moved.

UvA-DARE holds an open copy of the 1995 PDF at
`https://pure.uva.nl/ws/files/2908656/853_9632y.pdf`. It returned HTTP 503 on
four attempts across the session — an availability problem, not an access
problem. Worth one retry before asking the maintainer.

---

## 2. Correctly removed from `unclaimed` — 2 cases

### `J3.2` → `needs-paper`

The premise that its AIP DOI "almost certainly has nothing to do with it" is
wrong, and worth recording. `10.1063/1.4881599` is Pabst, *Analytical solution of
the Poisson–Nernst–Planck equations for an electrochemical system close to
electroneutrality*, J. Chem. Phys. **140**, 224113 (2014) — exactly this case's
structure (`S10`), and being *analytical* it would give a rank-2 validation under
the builder brief rather than a digitised figure. For a T0/P1 case that is worth
having.

But it is not the catalogue reference (Newman & Thomas-Alyea, a monograph not on
disk), and its only open copy is unreachable. The Jülich record is real —
JuSER record 154315, file `FZJ-2014-03663.pdf`, URL recovered from the OpenAIRE
index and recorded in the case file. `juser.fz-juelich.de` serves a JavaScript
bot challenge (a 248-byte `/fast-challenge/` stub) to curl and to WebFetch alike,
with or without a browser user-agent, cookie jar or referer. Working around that
is scraping past an access control, so it was not attempted; the page opens
normally in a browser. OpenAlex also flags that copy as `submittedVersion`, so
even once fetched it is an author manuscript whose equation numbering must be
checked against the published article.

Left `unclaimed`, this case would have sent a builder to a paper nobody holds.

### `J6.6` → `deferred`

`10.2139/ssrn.5388900` is Crossref `type: posted-content`, `subtype: preprint`,
`group-title: SSRN`, posted 2025, with no `relation` to a journal version and no
container title. AGENTS.md "Published work only" is unambiguous. `blocked_by` now
names the cheap automated test that lifts it: re-check Crossref for an
`is-preprint-of` relation on that DOI.

---

## 3. Metadata that turned out to be wrong

Four findings, all from cases nobody had opened.

| Case | Was | Actually | Action |
|---|---|---|---|
| `H1.9` | `10.1016/0009-2509(95)00102-b` | a packed-bed breakthrough paper with no membrane in it | DOI replaced with Krishna & Baur 2003 |
| `F1.6` | `10.1002/aic.690380821` | **a one-page book review**: Baird, "Bubble column reactors. By W.-D. Deckwer…", AIChE J **38**(8) 1305 (1992) | DOI removed |
| `F3.4` | `10.1149/1.2407312` | **a one-page book review**: "Gas-Liquid Reactions", J. Electrochem. Soc. **117**(10) 369C (1970), reviewer Arvo Lannus | DOI removed |
| `A4.7` | "Krishna & van den Broeke (1995)" | published author order is van den Broeke & Krishna | noted in the case file |

The two book reviews are a new failure mode and worth a rule: **when a case's
`catalog_reference` is a monograph, a resolved journal-article DOI is more likely
to be a review of that monograph than the monograph itself.** Both would have put
a reviewer's name on a gallery page in place of the author's. `J2.2`'s
book-chapter DOI was checked for the same reason and is *correct* —
`10.1016/b978-0-12-579650-7.50009-9` really is Randolph & Larson's chapter "The
steady-state MSMPR crystallizer", pp. 64–78.

---

## 4. The T0/P1 sweep — what it bought, and what it did not

66 cases sit at T0/P1 with `needs-paper`. Every one was put through Crossref,
Unpaywall, OpenAlex repository locations and — where a DOI looked Elsevier — the
authorised full-text API.

**Open-access routes found: three, none usable.**

- `J2.6` Smoluchowski (1917). The Jagiellonian Digital Library holds a free
  127-page scan. It was downloaded and opened: the cover sheet is Smoluchowski's
  own archive envelope, hand-lettered with the title, and p. 5 is the first
  manuscript page dated "Eingelaufen am 8. Sept. 1916" and annotated "Zeitschr.
  f. phys. Chem. 1917, 92, S 129–168". It is the **handwritten German
  manuscript**, not the printed article — 127 characters of extractable text in
  the whole file, with verso bleed-through on several leaves. Not kept.
- `J4.4` Luedeking & Piret. Biotechnology and Bioengineering republished it in
  2000 as a classic and both OpenAlex and Unpaywall flag that reprint as free
  (`publishedVersion`). Wiley's `pdfdirect` URL returned HTTP 403 to curl —
  Cloudflare, not a paywall notice. **This one should download from a browser in
  one click**; it is the cheapest item on the whole list.
- `J3.2` Pabst — covered above.

**Elsevier full text: 13 articles fetched, identity verified, zero cases
unblocked.** The reasoning is in `pdf-findings.md`; the short version is that
`content/object/pii` returns `{"choices":null}` for every one of them and the PDF
endpoint returns a one-page entitlement preview, so there is no page image to
render at 600 dpi — and the text itself loses decimal points on every scanned
original, up to and including 1991. Calderbank & Moo-Young returns 81 kB of text
containing five well-formed decimal numbers.

Promoting any of them to `unclaimed` on that basis would have been the exact
thing this repository exists not to do. They stay `needs-paper`, with three
things they did not have this morning: a **verified** citation, a **verified**
DOI, and a note saying which table a page would need — which is what makes the
request list below rankable instead of a guess.

**One genuinely new route did come out of it.** `A3.1` (Whitman two-film, 1923,
pre-DOI, previously unresolvable) is reprinted **verbatim** in Int. J. Heat Mass
Transfer **5**(5) 429–433 (1962) as "Pioneer Papers in Convective Mass Transfer,
5", with the header "W. G. WHITMAN: The two-film theory of gas absorption,
Chemical and Metallurgical Engineering 29, 146–148 (1923). Reprinted with
permission…". That reprint has its own DOI, `10.1016/0017-9310(62)90032-7`. This
is the AGENTS.md reprint route with an unusually clean provenance: cite Whitman
1923 as `reference`, the 1962 reprint as `reference_read_from`.

---

## 5. What to ask the maintainer for

Ranked. The maintainer presumably has one institutional login per publisher, so
the marginal cost of a second paper inside a batch is near zero — these are
ordered as **batches**, not as individual requests.

### Rank 1 — the ACS batch. 9 papers, 11 T0/P1 cases.

Nothing else comes close on cases-per-request. None of these has any API or OA
route.

| Paper | DOI | Unblocks |
|---|---|---|
| Froment, *Fixed bed catalytic reactors — current design status*, Ind. Eng. Chem. **59**(2) 18–27 (1967) | `10.1021/ie50686a006` | `C2.10`, `D3.4` |
| Dyson & Simon, *A kinetic expression with diffusion correction for ammonia synthesis on industrial catalyst*, IEC Fundam. **7**(4) 605–610 (1968) | `10.1021/i160028a013` | `C2.3`, `D3.2` |
| Fuller, Schettler & Giddings, *A new method for prediction of binary gas-phase diffusion coefficients*, Ind. Eng. Chem. **58**(5) 18–27 (1966) — **plus the erratum** | `10.1021/ie50677a007`, `10.1021/ie50680a601` | `A4.5` |
| Brunauer, Emmett & Teller, *Adsorption of gases in multimolecular layers*, JACS **60**(2) 309–319 (1938) | `10.1021/ja01269a023` | `J1.3` |
| Langmuir, *The adsorption of gases on plane surfaces of glass, mica and platinum*, JACS **40**(9) 1361–1403 (1918) | `10.1021/ja02242a004` | `J1.1` |
| Mears, *Tests for transport limitations in experimental catalytic reactors*, IEC Proc. Des. Dev. **10**(4) 541–547 (1971) | `10.1021/i260040a020` | `B1.7` |
| Voorhies, *Carbon formation in catalytic cracking*, Ind. Eng. Chem. **37**(4) 318–322 (1945) | `10.1021/ie50424a010` | `B2.1` |
| Chiu, Carratt & Soong, *A computer model for the gel effect in free-radical polymerization*, Macromolecules **16**(3) 348–357 (1983) | `10.1021/ma00237a002` | `J5.3` |
| Voltz, Morgan, Liederman & Jacob, *Kinetic study of carbon monoxide and propylene oxidation on platinum catalysts*, IEC Prod. Res. Dev. **12**(4) 294–301 (1973) | `10.1021/i360048a006` | `C2.13` — **but try the reprint route first**: Oh & Cavendish (1982) is already on disk as the `I1.2` file. Ask the `E1.1` question of it — does it *name and test* the Voltz expression, or merely use it? |

ACS scans have historically had the best text layers of anything in this
repository (`i260028a001.pdf` at 12.4 k chars/page, `i300005a006.pdf` at 9.0 k),
so this batch is also the cheapest to extract once it arrives.

### Rank 2 — the Wiley / AIChE Journal batch. 7 papers, 7 T0/P1 cases.

Wiley has no API route at all, so every one of these is blocked on a manual pull.

| Paper | DOI | Unblocks |
|---|---|---|
| Bischoff, *Effectiveness factors for general reaction rate forms*, AIChE J **11**(2) 351–355 (1965) | `10.1002/aic.690110229` | `B1.3` |
| Myers & Prausnitz, *Thermodynamics of mixed-gas adsorption*, AIChE J **11**(1) 121–127 (1965) | `10.1002/aic.690110125` | `J1.4` |
| Heck, Wei & Katzer, *Mathematical modeling of monolithic catalysts*, AIChE J **22**(3) 477–484 (1976) | `10.1002/aic.690220310` | `I1.3` |
| Larkins, White & Jeffrey, *Two-phase concurrent flow in packed beds*, AIChE J **7**(2) 231–239 (1961) | `10.1002/aic.690070213` | `G1.1` |
| Bhatia & Perlmutter, *A random pore model for fluid–solid reactions: I*, AIChE J **26**(3) 379–386 (1980) | `10.1002/aic.690260308` | `B3.3` |
| Andrews, *A mathematical model for the continuous culture of microorganisms utilizing inhibitory substrates*, Biotechnol. Bioeng. **10**(6) 707–723 (1968) | `10.1002/bit.260100602` | `J4.2` |
| Yagi & Kunii, *Studies on effective thermal conductivities in packed beds*, AIChE J **3**(3) 373–381 (1957) | `10.1002/aic.690030317` | `A3.12` — **check the reprint route first**: the same authors' 1955 CES paper is on disk as the `B3.1` file. |

**Free bonus in the same session:** Luedeking & Piret's 2000 Biotechnol. Bioeng.
reprint, `10.1002/(SICI)1097-0290(20000320)67:6<636::AID-BIT3>3.0.CO;2-U`, is
flagged open access and downloads from a browser. Unblocks `J4.4`. One click.

### Rank 3 — Elsevier PDFs. 8 papers, 7 T0/P1 cases, lowest risk of a wrong pull.

Each of these has already been fetched as text, so its identity is confirmed and
the exact table a page needs is known. Only the page images are missing.

| Paper | DOI | Unblocks | Why it is worth it |
|---|---|---|---|
| Uppal, Ray & Poore, *On the dynamic behavior of continuous stirred tank reactors*, CES **29**(4) 967–985 (1974), **plus** the 1976 sequel CES **31**(3) 205–214 | `10.1016/0009-2509(74)80089-8`, `10.1016/0009-2509(76)85058-0` | `D2.4` | **Best page of the lot.** The model is the standard dimensionless non-adiabatic CSTR and the validation is rank 2 in the builder brief — the paper's own analytic multiplicity and limit-cycle criteria, independently re-derived and checked against numerical continuation. A check that can fail. The numerical cases are already half-legible (`B=15.0 BETA= 0.40`, `Da=0.0006017` survive in the figure labels) so the PDF mostly buys the criteria. |
| Mars & van Krevelen, *Oxidations carried out by means of vanadium oxide catalysts*, CES **3** (Spec. Suppl.) 41–59 (1954) | `10.1016/S0009-2509(54)80005-4` | `C1.3` | Five tables of kinetics — rank 1/3 validation, no figure digitisation. The text version is unreadable (four different glyphs for one decimal point). |
| Hulburt & Katz, *Some problems in particle technology: a statistical mechanical formulation*, CES **19**(8) 555–574 (1964) | `10.1016/0009-2509(64)85047-8` | `J2.1`, and underwrites `J2.2` and `J5.1` | The origin of the population balance and its method of moments. Broadest reach of any Elsevier item. |
| van Deemter, Zuiderweg & Klinkenberg, *Longitudinal diffusion and resistance to mass transfer as causes of nonideality in chromatography*, CES **5**(6) 271–289 (1956) | `10.1016/0009-2509(56)80003-1` | `J1.10` | Two experimental data tables, both destroyed in the text (`0 163 1100 0 2? 0 30`). The API bought this case nothing. |
| Robeson, *The upper bound revisited*, J. Membr. Sci. **320**(1–2) 390–400 (2008) | `10.1016/j.memsci.2008.04.030` | `H1.8` | **Ask for the 2008 paper, not the 1991 one.** It is born-digital and will extract cleanly; the 1991 original is a scan whose Table 2 loses every decimal point, and the API returns 404 for the 2008 PII. |
| Szekely & Evans, *A structural model for gas–solid reactions with a moving boundary — II*, CES **26**(11) 1901–1913 (1971) | `10.1016/0009-2509(71)86033-5` | `B3.2` | Part II carries the porosity/grain-size/temperature parameter study; Part I is mostly formulation. Ask for II, or both. |
| Calderbank & Moo-Young, *The continuous phase heat and mass-transfer properties of dispersions*, CES **16**(1–2) 39–54 (1961) | `10.1016/0009-2509(61)87005-X` | `A3.6` | 81 kB of text, five readable numbers. |

### Rank 4 — monographs

Restating the standing item from `handoff.md`, now with the specific cases this
sweep tied to each:

- **Newman & Thomas-Alyea, *Electrochemical Systems*** → `J3.2` (T0/P1, `S10`).
  Or, as an alternative that gives a stronger validation, the published PDF of
  Pabst, J. Chem. Phys. **140**, 224113 (2014), `10.1063/1.4881599`.
- **Randolph & Larson, *Theory of Particulate Processes*** (1971) → `J2.2`; the
  DOI `10.1016/b978-0-12-579650-7.50009-9` is confirmed to be its chapter 4, and
  Elsevier's API does not serve book chapters.
- **Deckwer, *Bubble Column Reactors*** → `F1.6`, which now has no DOI at all.
- **Danckwerts, *Gas–Liquid Reactions*** (1970) → `F3.4`, likewise.
- Bird/Stewart/Lightfoot, Taylor & Krishna, Levenspiel, Froment & Bischoff — the
  ~20-case bundle already on the dashboard.

### Not worth asking for

- **`J2.6` Smoluchowski.** The free copy is the handwritten German manuscript;
  the printed article is at De Gruyter behind a subscription. Low value per unit
  of effort even if obtained.
- **`D1.1`–`D1.5`.** These five T0/P1 cases have a `catalog_reference` that is a
  *structure code* (`S2`, `S4`, `S7`, `S8`, `S6`), not a citation. No paper will
  unblock them; they need a scope decision about which exemplar paper each should
  reproduce. They should not be sitting in `needs-paper` at all.

---

## Appendix — what was actually downloaded

One PDF kept (`A1.8`, MFIX Theory Guide, OSTI). One PDF fetched, verified and
discarded (`J2.6`, Smoluchowski manuscript, 94 MB, handwritten). 15 Elsevier
full texts fetched, all identity-verified and kept under
`~/papers/pymrm-gallery/api-text/` with their per-file quality notes in
`pdf-findings.md` and in each case file. Nothing was fetched from behind a
paywall, a login or a bot challenge; the three routes that ended in one
(Wiley `pdfdirect`, JuSER, UvA-DARE) are recorded with their exact URLs so a
human can finish them in a browser.
