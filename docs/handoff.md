# Handoff — state as of 2026-08-03

Start here if you are picking this up fresh. Read this, then
[`AGENTS.md`](../AGENTS.md), then [`pdf-findings.md`](pdf-findings.md).

**Read [“The check that cannot fail”](#the-check-that-cannot-fail--the-defect-the-verifier-exists-to-catch)
before writing any validation cell.** Across 2026-07-31 to 08-03, twenty-four
pages went through adversarial verification and **eighteen** carried some form of
it — a check that could not fail, a claimed sensitivity the check did not have,
or an independence claim that did not survive checking. It is by a wide margin
the most common finding in this repository, and it has now been found **five
times inside the break table or diagnostic built to guard against it**: `A3.4`,
`D2.1`, `A1.5` (14 of 36 rows in the headline wall-law fit were an identity),
`A1.8` (every check sat on the dense branch, so two constants could be deleted
outright), and `A4.7` (the sole diagnostic differentiated the same object it
tested, so it was guaranteed for any parameter values).

**A second class was found and swept on 2026-08-02**: a page that loads another
page's CSV and does not read that page. See the rule in
[`AGENTS.md`](../AGENTS.md) under Data rules, and
[`cross-page-audit-2026-08-02.md`](cross-page-audit-2026-08-02.md). It reached a
reported number on `J3.1` and the Reuse advice on `A1.6`.

## Built and live

**https://computational-chemical-engineering.github.io/pymrm-gallery/**

40 page directories, 41 published catalog entries, both CI workflows green.
(The two counts differ because `B1.1` covers `B1.5` and `A2.1` covers `A2.2`.)

### Added 2026-08-02 to 08-03 — six pages

| Page | What it shows | Validation |
|---|---|---|
| `A1.5` Richardson–Zaki | `u/u_t = ε^n` is not what the paper writes; the wall term is the result | 66 slopes off 600 dpi renders, 5 printed errors found; eq. (34) **is** its window's mean, 2.390000 vs printed 2.39 |
| `A1.6` Wen & Yu u_mf | Where 33.7 and 0.0408 come from, and what they cost | Held out against Geldart's 21 velocities: −17.4 %; at *his* voidage the reference balance is +58 % — worse than what it judges |
| `A2.5` Edwards–Richardson | How much of a fitted dispersion coefficient is the grid | **experimental** — 55 points digitised by computing and *erasing* the printed curve; null baseline 30.3 % |
| `A1.8` Gas–solid drag closures | Two closures agreeing in both limits and not in between | Factor **1.73** at Re_m = 4.09 where the limits are 1.192 and 0.957 |
| `A4.7` Zeolite micropore M–S | The M–S form removes a *divergence*, not a dependence | 8 printed targets, worst 0.78 %; the factor-of-two confirmed out-of-sample on CH₄ |
| `B1.2` Aris (repair) | No rescaled modulus explains Table 1's two cells | 0.197 needs Λ×0.961, 0.100 needs ×0.974 — two slips, and the page says it is an inference |

### Added 2026-07-31 to 08-02 — twenty pages

| Page | What it shows | Validation |
|---|---|---|
| `A4.2` Maxwell–Stefan vs Fick | Two bookkeepings of one physics; the scalar closure fails *structurally* | MS invariant to species ordering to 3.9e-16 where Wilke moves 1.44–4.31 mole % |
| `J3.5` SPM / SPMe | Asymptotic reductions of `J3.4`, and what their slope test resolves | Tables II+III 21/21; the exponent drifts 2.07→2.70, stated as a window |
| `H1.4` Itoh Pd membrane reactor | Conversion past equilibrium — and how little the measurement tests | Stated 99.7 %; **45 % low on 1−X**, the quantity that varies |
| `F1.3` Wilkinson holdup | A correlation whose one pressure term fails | 6.8 decades low at helium — the authors' own published finding |
| `E1.2` Davidson bubble | Eqs. 9 and 10 are one flow field read twice | The printed 4.5 recovered to 1.1e-6, never an input |
| `B1.6` Prater relation | An identity, and the two ways it breaks | Defect-sensitivity table; film and Le ≠ 1 breakdowns |
| `J3.1` Butler–Volmer | The law, and where each limit stops being safe | Eq. 17 collapses to the separately-printed Eq. 30, factor for factor |
| `B1.4` Weisz–Prater criterion | The inverse problem: what an observed rate can tell you | Non-uniqueness *inside* the safe band — Φ = 0.32 admits η = 226 |
| `H1.1` Sieverts permeation | The half-power law in its own right | α_H derived, not quoted; Figure 1's printed bore breaks the circularity |
| `A4.4` Knudsen / Bosanquet | Bosanquet is exact under Knudsen control, worst under bulk | `D_Bos/D₁ = 1 + (ρ−1)x₁/(1+Kn)`, re-derived independently |
| `D2.1` Barkelew diagram | The collapse drifts | τ_m 1.2879→1.1622 over S = 4–200, Richardson-extrapolated |
| `A4.3` Dusty gas model | The n-component form and the viscous term | Wall term proved to make [B] invertible; Haynes's remark holds everywhere |
| `A1.1` Ergun | Four correlations against the data Ergun fitted | **experimental** — k₁ = 151.9, k₂ = 1.697 vs his 150 and 1.75 |
| `G1.7` Wammes–Westerterp | High-pressure trickle-bed hydrodynamics | **experimental** — 2.9 % holdup vs the paper's own 8 % |
| `F3.5` CO₂–amine absorption | A reconstruction, with the evidence against it reported too | J_u within 3.8 %; desorption *favours* the printed constants |
| `G1.8` Trickle-bed partial wetting | The legend is offset one row against the curves | Separability forces the spacing; three curves land on printed values |
| `A2.1` Danckwerts BCs (covers `A2.2`) | The canonical outward-normal case | Closed form at order 2.00; the outlet trap is *undetectable* |
| `A1.7` Geldart classification | The boundaries, and what his own table says | **experimental** — but a null predictor scores 19/21, so 2-of-3 is the claim |
| `A3.4` Wakao–Funazkri | The correlation, and what the figure can and cannot settle | α = 1.100 with β fixed; the free fit is estimator-dependent |
| `B1.2` Aris shape modulus | The general case Aris called "excessively difficult" | 14 shapes; the sphere is **not** the floor — the cube lies 3.7 % below |

Also corrected on already-published pages: `F1.4`'s mis-quoted measured band,
`A4.2`'s provenance overclaim, and a dense-Jacobian trap on `B1.1`, `B1.6` and
`F3.1` (6.3× and 2.8× runtime, answers bit-identical).

Eight of these came from the **reprint route** — the origin paper unreachable, a
paper on disk printing it in full with attribution. See `reference_read_from` and
`origin_not_consulted` in [`AGENTS.md`](../AGENTS.md).

| Page | What it shows | Validation |
|---|---|---|
| `A4.9` Duncan–Toor ternary diffusion | Osmotic → uphill → diffusion barrier | **experimental** — 0.59 mole % vs the paper's own 0.45%, 28 digitised points |
| `C2.1` Xu–Froment steam reforming | The most-used SMR kinetics, against the runs they were fitted to | **experimental** — 0.0017 in conversion (2.7%) over 61 digitised points, nothing fitted |
| `B1.1`+`B1.5` Thiele + Weisz–Hicks | η(φ), and 3 steady states at one φ | 2.2e-4 vs exact; both methods agree on all three branches |
| `D2.2` Van Welsenaere–Froment runaway | Two criteria for the runaway boundary, swept over the operating plane | 0.054% over all 30 numbers in their Section 6; two independent methods agree to 0.18% |
| `F3.1` Hatta regimes | Enhancement factor, 3 regimes | 6.3e-3 vs exact; VKH good to 2.1%, DeCoursey to 8.7% |
| `A2.3` Taylor–Aris dispersion | Homogenisation closure, and when the lumped coefficient becomes defensible | 1.0e-4 vs Eq. 25 at n_r=200, O(h²); Taylor's own capillary run to 0.04% |
| `J1.5` LDF breakthrough | What the linear-driving-force constant actually is | 6.6e-5 vs the exact series at n_r=400 |
| `F1.4` Krishna–Ellenberger holdup | A correlation with no fluid property in it, against the figure it was fitted to | **experimental** — 13.8% mean deviation, +2.8% bias over 63 digitised points |
| `H1.7` Wijmans–Baker solution–diffusion | Two constants predict the third — and the figure cannot test the prediction | **experimental** — A and B fitted to 8 points, rejection predicted for the other 4 |
| `B3.1` Yagi–Kunii shrinking core | The one equation all three textbook regimes come from, plus a map of where they are safe | 6.9e-16 against an independent derivation; all three limits to 2.4e-8 |
| `F2.3` Maretto–Krishna FT slurry column | Plug-flow large bubbles over a well-mixed slurry, and two printed constants that stop it working | **experimental** holdup — 5–6% over 79 points; conversions 93.1/63.8% vs the paper's 96/63% |

Status counts live in `models.yaml`: 22 published, 12 planned, 2 deferred
(`H1.12`, `B1.12` — unpublished manuscripts, see the published-work-only policy
in [`blueprint.md §9`](blueprint.md#published-work-only-policy)).

**Four pages are validated against experiment** (`A4.9`, `C2.1`, `F1.4`, and
`F1.3` on `F1.4`'s dataset), plus `A2.3` against Taylor's own worked capillary
run and `H1.4` against a stated measurement — though `H1.4`'s page now says
plainly how little that comparison resolves. Everything else is tier 6 by
necessity: the source papers contain no measurements. `A3.4` is still the next
chance to move the ratio and still needs figure digitisation.

### The queue, and why it is paper-starved

266 catalogued cases as of 2026-08-03: **41 published**, 5 covered, 3 deferred,
1 unclaimed, and **216 blocked on a source PDF**. Nothing is waiting on a
maintainer decision — `standing-decisions.yaml` is empty and every parked figure
review has been answered.

That 216 is irreducible by automation, and this has now been established twice
rather than assumed. `find_papers.py` over all 266 found five open-access routes;
a targeted sweep of the 66 T0/P1 cases on 2026-08-02, with the authorised
Elsevier key, unblocked **zero** — see *Papers available* for why. Most of
section A predates DOIs entirely.

So the productive veins are the **reprint route**, cases whose origin turns out
to need no paper at all (`A1.8`), and — the largest single lever, still unspent —
**one monograph on disk** for the approved textbook-canonical class.

## Papers available

**Thirty PDFs** at `~/papers/pymrm-gallery/`, and **every one has now been
consumed by a published or parked page.** Inventory with per-file text-layer
quality is in [`pdf-findings.md`](pdf-findings.md); the catalogue-ID → filename
map is [`papers-on-disk.yaml`](papers-on-disk.yaml) and must be updated by hand
whenever a PDF arrives — automatic matching misses publisher-PII filenames, and a
mis-mapped entry asks the maintainer for a paper they already supplied (`H1.4`
was listed under `H1.9` until 2026-07-31).

**A filename is not a paper — open it before you believe it.** This has now cost
four separate errors. Two auto-resolved DOIs point at *book reviews* of the
monograph the catalogue names (`F1.6` → a one-page review of Deckwer, `F3.4` →
one of Danckwerts); `H1.9`'s recorded 1995 PII is a packed-bed carbon paper with
no membrane in it; and on 2026-08-02 an agent read
`1-s2.0-0009250957850283-main.pdf` as a scan of Aris 1957 and dispatched work on
that basis. **It is one page** — the Elsevier entitlement preview, journal page
262 — and `B1.2`'s open question is on page 265. `papers-on-disk.yaml` now
carries a loud comment saying so. New rule worth keeping: when
`catalog_reference` is a monograph, an auto-resolved *journal* DOI is likelier to
be a review of the book than the work itself.

**Do not expect the Elsevier API to unblock a pre-1995 case.** Swept across 66
T0/P1 cases on 2026-08-02 and it unblocked **zero**. The key reaches the text,
but `content/object/pii` returns no page images and the PDF endpoint returns a
one-page preview, so **there is nothing to render at 600 dpi** — and the text
drops decimal points on every scan through 1991 (Robeson's Table 2 comes back as
`2 6 2 89 3 3 3 46 3 64 3 8`). Correcting an earlier caveat: **1995 is not
"post-scan era" either** — that PII corrupts species labels and omits its Table 4
entirely. Krishna & Baur 2003 *is* born-digital and clean, which is why it
sources `A4.7`. The ranked, publisher-batched request list is in
[`source-sweep-2026-08-02.md`](source-sweep-2026-08-02.md).

Three arrived during the 2026-07-31 sweep: `Marquis_2019…` (`J3.5`, gold OA from
IOP), `F3_5.pdf` and `G1_7.pdf` (both from the University of Twente repository).
Note both of the latter were **auto-resolved to the wrong citation** and the
builders' identity checks caught it — `F3.5` is Bosch et al. 1989, not the
catalogue's Versteeg & van Swaaij 1988; `G1.7` is the Chem. Eng. Technol. paper,
not the CES one. Always verify identity before building.

## Elsevier full text — use this instead of OCR

An institutional API key is stored at **`~/.config/elsevier/apikey`** (mode 600,
outside the repo; `.gitignore` blocks credential filenames — never commit it).

The **DOI endpoint 404s** on older Chemical Engineering Science articles because
their DOIs contain parentheses. **Use the PII endpoint**, and note the PII is
recoverable straight from the ScienceDirect filename:
`1-s2.0-`**`0009250970850734`**`-main.pdf`.

```bash
K=$(cat ~/.config/elsevier/apikey)
curl -sS -H "X-ELS-APIKey: $K" -H "Accept: text/plain" \
  "https://api.elsevier.com/content/article/pii/0009250970850734"
```

Verified working on Van Welsenaere & Froment (47 kB) and Wakao & Funazkri
(45 kB) — both of which had the *worst* local OCR of the set.

**Correction, 2026-07-26 — do not trust the API text for numbers.** An earlier
version of this file claimed the API "returns properly encoded text, so
exponents and subscripts survive". That is wrong for pre-1980 scans. What comes
back is the *publisher's* OCR of the same scan, and on Van Welsenaere & Froment
it drops decimal points wholesale:

| API text | Truth (600 dpi page image) |
|---|---|
| `R = 00125 m` | *R* = 0.0125 m |
| `M = 2948 kg/kmole` | *M* = 29.48 kg/kmole |
| `(p°),,, = 001353 atm` | (*p*⁰)ₗ = 0.01353 atm |
| `b = 19837` | *b* = 19.837 |
| `t,„ = 21 -818` | *t*ₘ = 21.818 |

The 1970 typesetting uses a mid-dot decimal separator, which the OCR discards.
So the API is excellent for **prose** — section structure, what each figure
shows, which assumptions are stated — and useless for **parameters**. Read
those from a 600 dpi page render, the same as for a Wiley scan. Treat any
integer with an implausible magnitude as a lost decimal point and go to the
image; never "fix" it by inference.

**Other limits:** Elsevier only. Xu & Froment and Krishna & Ellenberger are
Wiley and have no API route at all. And full text never yields figure data —
plots still need digitising. Be polite with request volume; this is not a
bulk-download tool.

## Extraction cost, revised

| Page | Paper | Extraction | Value |
|---|---|---|---|
| `D2.2` | Van Welsenaere & Froment | API for prose, **page image for every number** | Runaway criteria; sweep-based figures are striking |
| `A3.4` | Wakao & Funazkri | API for prose, page image for numbers | Sh–Re dataset, but likely a scatter *figure* → digitise |
| `E2.1` | Kunii & Levenspiel | good OCR (12.4k chars/page) | The canonical fluidised-bed model |
| `I1.2` | Oh & Cavendish | good OCR (9.0k) | Converter light-off, `S4`+`S7` |
| ~~`F1.4`~~ | ~~Krishna & Ellenberger~~ | ~~tables clean, holdup in figures~~ | **done 2026-07-28** |
| ~~`C2.1`~~ | ~~Xu & Froment~~ | ~~hard~~ | **done 2026-07-26** |
| ~~`D2.2`~~ | ~~Van Welsenaere & Froment~~ | ~~API + page image~~ | **done 2026-07-27** |

## Recommended next moves

**The section-A batch is finished and the queue is now genuinely paper-bound.**
216 of 266 cases need a PDF, 41 are published, 5 are covered, 3 deferred, and at
the time of writing exactly **one** case (`H1.9`) is buildable from material in
hand. Every PDF on disk has been consumed.

1. **Ask for papers in publisher batches, not one at a time.** The marginal cost
   of a second paper inside one login is near zero, and
   [`source-sweep-2026-08-02.md`](source-sweep-2026-08-02.md) ranks them that
   way. The ACS batch is the best single ask: **9 papers unblocking 11 T0/P1
   cases**, and ACS scans have the best text layers in the repo. One item is
   free — Luedeking & Piret's 2000 Biotechnol. Bioeng. reprint is open access and
   downloads from a browser (curl gets a Cloudflare 403), unblocking `J4.4`.
2. **The textbook-canonical approval is banked but inert.** The ~20 T0 cases are
   approved to build without the original paper, but `AGENTS.md` still forbids
   writing from memory, and **none of the monographs is on disk**
   (Bird/Stewart/Lightfoot, Taylor & Krishna, Levenspiel, Froment & Bischoff).
   One book unlocks more than any paper on the dashboard. Until then the class
   cannot move.
3. **Five cases do not belong in `needs-paper` at all.** `D1.1`–`D1.5` are T0/P1
   with a *structure code* rather than a citation as their `catalog_reference`.
   No paper unblocks them; they need a scope decision.
4. **Work the reprint route.** Nine pages came from it. The test, from
   `AGENTS.md`: a paper on disk that prints the result *and* names or tests it.
   Two verdicts bound it — `E1.1` failed because Kunii & Levenspiel print the
   two-phase relation but never name, attribute or test it; `J3.3` *passed* the
   test and was still `covered`, because the theory is J3.4's structure. Ask both
   questions, in that order. One new route is already located and unused:
   **`A3.1` Whitman 1923 is reprinted verbatim in IJHMT 5 (1962) 429–433** under
   its own DOI.
5. **Some cases need no paper ever.** `A1.8`'s origin is an *unpublished* 1987
   report, so the DOE Theory Guide on disk is the citable published source and
   the case is closed rather than parked. Check the reference list of what you
   have before adding to the request pile.

### When a printed constant is wrong, prove it from the paper's own results

`F2.3` needed two corrections before it would run, and the method for
establishing them generalises.

Eq. 2's rate prefactor is printed as 8.8533e3 mol/(s kg_cat bar²), which gives an
intrinsic rate 10⁶ larger than any cobalt catalyst. Eq. 1's rate is labelled
`R_CO+H2` but behaves as a CO rate. Neither could be fixed by fitting — that would
have made the whole comparison circular.

**What made the diagnosis safe was establishing chemical control first.** The
paper reports that a 10-fold rise or 3-fold fall in kLa is negligible; reproducing
that *before* touching the kinetics proves the mass-transfer correlations are not
free to absorb a rate error. Only then does the conversion comparison isolate the
kinetics, so each correction becomes a discrete choice between stated
alternatives — 10³ vs 10⁻³, syngas vs CO — with the paper's own reported
conversions selecting between them. The page prints the alternatives and what each
gives.

*Order matters: pin down what is NOT free before claiming a constant is wrong.*

Two pymrm traps recorded on that page, both silent failures:

- A convection outlet left as `None` makes the matrix **singular**, and a
  rank-deficient solve still returns a plausible-looking profile.
- With varying velocity, discretise `d(Uc)/dz` as the divergence of the flux.
  `U dc/dz` loses the gas contraction — 65 % of the volumetric flow here.

### The check that cannot fail — the defect the verifier exists to catch

**Read this before writing a validation cell.** On 2026-07-31 six pages went
through adversarial verification and *four* carried the same defect, in different
disguises. It is now the most common finding in this repository, it always looks
like the page's strongest evidence, and inspection never catches it.

The shape: two routes to a quantity are compared, they agree to machine
precision, and the agreement is **algebraically guaranteed**, so the check has no
power against the error class it is presented as guarding.

- `A4.2` compared Maxwell–Stefan against generalized Fick at **8.9e-16**. Both
  called the same `build_b`; one solved the 2×2, the other wrote out its
  adjugate. A dropped sign inside `build_b` moved both together and the metric
  still read machine precision. Replaced by the n×n friction system, which never
  forms [B] — that one agrees to 7.8e-16 and *can* fail.
- `F3.5` reported electroneutrality at 2.8e-11 and charge flux at 1.1e-9. All
  charged species share one diffusivity, the reactions conserve charge and both
  boundaries are Dirichlet: those residuals cannot be anything else.
- `F1.3` "recovered" a printed exponent 0.757 as 0.7566 — from a function
  containing `**0.757`. Substituting 0.657 recovers 0.6557. It was claimed in
  five places as the defence against a mis-read exponent.
- `J3.5` claimed a slope test resolved four correction terms; deleting one of
  them entirely changed the answer by 0.02 %.

**`B1.6` is the worked example, because there the break test was actually run.**
Its identity residual ε held at ~1e-11 across 324 solves. Injecting defects showed
what that number is blind to:

| injected defect | result | ε |
|---|---|---|
| `maxfev=1`, Newton residual 18 | centre y = **2.12** (impossible, y > 1) | 8.5e-12 |
| sphere solved with `nu = 0` | η **57 % wrong** | 4.7e-12 |
| `n_u = 3` | η **37 % wrong** | 2.4e-15 (*better* than at 200) |
| garbage 5-node mesh in the "independent" `solve_bvp` route | — | 2.0e-11 |

So ε could not see an unconverged solve, a wrong geometry, a wrong grid, or a
garbage mesh — and the page claimed all four. The geometry claim was impossible by
construction: one `construct_div` serves both fields, so geometry cancels out of
the invariant identically. The independent discretisation was not independent
*for this quantity*: `w = θ + βy` is a closed linear subsystem that scipy
inherits exactly.

What ε *does* catch, measured the same way: a source-sign error → 0.42, a 1 %
source-scale error → exactly 1.0e-02, a boundary-condition-type mismatch → 0.16,
an inconsistent rate state → 0.54. That is a real and useful check — it detects an
inconsistency between the two source terms or the two boundary conditions, and
nothing else. **The fix was never to delete the check; it was to say what it
tests.** Publishing the sensitivity table is strictly better than publishing the
residual alone.

**A powerless check hides a second failure: nobody notices the reference is
broken.** Fixing `B1.6` turned up something the verifier had not reached — **6 of
the 12 `solve_bvp` collocation runs never reached their own tolerance**
(`status = 1`, rms residual up to 4.0e-01), and the reported "worst 2.6e-13" came
from one of the failed ones. The invariant was satisfied by a solution that was
not a solution. Adding the check that *can* fail — pymrm profiles against
collocation profiles — both converged at order 2.00 on the good references and
stalled at 9.3e-05 with order −0.00 on the bad one, flagging it immediately.
**Assert your reference solver converged; do not infer it from an identity.**

**Two pages carried a code comment that asserted a sensitivity the check did not
have** (`B1.6`, `E1.2` — "a wrong nu or a mis-signed boundary flux breaks it
immediately"; measured, it moves for neither). The comment gets written when the
check is conceived, the check weakens as the code evolves, and nobody re-runs the
claim. **A comment claiming sensitivity is a claim, and needs the same break test
as the number.**

Three questions to ask of every agreement number, before it goes on a page:

1. **Do the two routes share code?** Same assembly, same operator, same
   parameter dict — then you are testing arithmetic, not physics.
2. **Deliberately break something the check should catch** — wrong `nu`, a
   flipped sign, a mismatched boundary condition — and confirm the number
   *moves*. If it does not, the check is decoration. This is cheap and it is the
   only reliable test.
3. **Is the identity structural?** Conservation residuals are often exact by
   construction. Say what the check confirms (bookkeeping, an implementation
   port) rather than implying it confirms the model.

A check that cannot fail is not evidence, and presenting one as evidence is the
kind of quiet overclaim that costs more credibility than a missing page. Keeping
it is fine — label it for what it is, as `A4.2` now does ("algebraic identity …
cannot detect an error inside [B]").

**The sibling defect: prose numbers drift from code output.** Four pages in the
same batch printed a number in markdown that the notebook contradicted two cells
above — 30 % against 104.8 %, "twenty" against 21, "15–20 %" against 8.4–21.4 %,
"6 to 14 Newton iterations" against 6–27. Interpolate computed values into the
prose instead of typing them, and before reporting `ready`, re-read every
markdown number against what the cells actually print.

### The reprint route: what it has produced, and where it stops

Every paper on disk is consumed, and 226 cases have no reachable source — so the
only route still open is a paper on disk that prints someone *else's* result.
`AGENTS.md` states the test; this is what running it eight times taught.

**It works, and it is not marginal.** `B1.6` (Prater, via Weisz & Hicks),
`F1.3` (Wilkinson, via Krishna & Ellenberger), `E1.2` (the Davidson bubble, via
Kunii & Levenspiel), `J3.1` (Butler–Volmer, via Doyle), `B1.4` (Weisz–Prater,
via Weisz & Hicks), `H1.1` (the half-power law, via Itoh) and `A4.4` (Bosanquet,
via Krishna & Wesselingh) all came from it.

**Two verdicts bound it, and both are worth knowing before dispatching.**

- `E1.1` → **needs-paper**. Kunii & Levenspiel print the two-phase relation but
  never name, attribute or test it. A reprint that merely *contains* the equation
  cannot source a page *about* the result.
- `J3.3` → **covered**. Doyle names, attributes and tests porous electrode
  theory — the test passes comfortably — but the theory *is* J3.4's structure,
  and the two places it is open to examination are where J3.4 already spends its
  validation. Passing the reprint test does not entitle a case to a page.

So ask two questions, in order: *does the source name and use the result?* and
*does an existing page already examine it?* The second killed a case the first
had cleared.

**Check the source carries the case before dispatching, not after.** A `grep`
costs seconds. `J4.1` (Monod) was dropped from a batch because ASM1 mentions
Monod exactly once — the E1.1 failure caught at the cheapest possible moment. I
also dispatched `E1.1` on a premise that turned out false, and its builder
refuted me by reading the citation list; and briefly claimed `F3.4` before
noticing its source is a book that is not on disk.

**Expect the received attribution to be wrong.** Three of the last four cases
found the printed page contradicting the name the catalogue uses:

| case | catalogue says | the source actually says |
|---|---|---|
| `B1.4` | Weisz & Prater 1954 | "It was shown by WEISZ **[9]**" — Z. Phys. Chem. 1957 |
| `E1.2` | Davidson & Harrison | "the Davidson bubble", once, in a list of assumptions; no entry in the literature cited |
| `H1.1` | Sieverts | half-power law attributed to Bohmholdt & Wicke 1967; Sieverts & Danz cited only for C₀ |

None of these changes what the page builds — eq. (1) really is the Weisz–Prater
group — but each changes what the page may *claim*. Record the origin under
`origin_not_consulted:` (the `J3.1` form) when the source cites rather than
reprints it, and `reference_read_from:` when it genuinely reprints the result.

### The Neumann outflow extrapolates to the face — a hand-written outlet flux misses

Found on `A3.7` (2026-08-05), confirmed by its verifier and reproduced again by
its fixer. pymrm's zero-gradient outflow boundary **extrapolates the value to
the face** rather than taking the upwind cell, so a hand-written outlet flux of
the form `v·C_N` (last cell centre) disagrees with what the operator actually
transports: at n = 8 the mismatch is **8.9e-3**, falling first-order to 7.8e-5
at n = 800. A mass balance written with `v·C_N` therefore fails to close by
~1e-4 on a modest grid **and looks exactly like a physics error**. Read outlet
values through `compute_boundary_values` (or the face itself), never off the
last cell centre — the same lesson as `A2.6`'s outlet metric, which was
published 11.4 % low for reading h/2 short of the boundary.

### `NumJac(shape)` on a one-field 1-D problem builds a dense Jacobian

Found while building `B1.4`, confirmed independently, and it was live on three
published pages. The default stencil couples the **last** axis in full, so a bare
`(n,)` shape declares every cell coupled to every other:

| n | `NumJac((n,))` constructor | nnz/n² | `NumJac((n,1))` |
|---|---|---|---|
| 400 | 3.2 s | 1.0000 | 0.0005 s |
| 1600 | 70 s | 1.0000 | 0.0012 s |

Whole-page effect, answers **bit-identical** in every metric: `B1.1` 321 s → 52 s
(6.2×), `B1.6` 64 s → 22 s (2.8×). Always write the single-field layout as
`(n, 1)` so the last axis is the field axis — which is the house convention
anyway, and degrades gracefully when a second field appears.

**And do not reach for `axes_diagonals=[0]` on a 1-D shape.** `AGENTS.md` used to
say to add it "when the source term depends on neighbouring cells", with no
dimension caveat. On `ndims=1` that is *wrong*, not just wasteful: `axes_blocks`
still defaults to `[-1]`, the same axis 0, so the `[-1, 0, 1]` offsets are
reinterpreted as absolute indices n−1, 0, 1. The Jacobian comes out with no
diagonal and the solve converges to a different answer. It is only meaningful at
`ndims ≥ 2` — `NumJac((n, 1), axes_diagonals=[0])` — and only when the *source
term itself* reads neighbours; the Laplacian's coupling normally arrives
analytically through the divergence operator.

*The general lesson: a convention documented without its precondition is worse
than no convention, because agents follow it. When a rule turns out to be
dimension-dependent, say so where the rule is written.*

### A curve the tracker never saw is invisible to every check built on it

`F3.5`'s maintainer review returned two corrections — Figure 2 was missing its
top solid curve, and Figure 3's two carbonates had been fused into one. Both had
the same root cause, and the useful part is *why the extraction could not catch
it itself*.

The missing Fig. 2 curve **leaves the plot box through the top**, ending above
the highest tick, so a tracker whose read band starts at the frame never sees it.
The Fig. 3 pair converge to within one line width at the frame — a single 8 px run
where one rule is 6–7 px — so they merged. Neither failure produces a *wrong*
number; both produce a *missing* one, and nothing in the extraction complains.

**The check that would have caught it could not be written until the curves were
right.** Charge balance on the promoted solution needs the promoted carbonate of
every figure — exactly the curve that was missing. Once added, all three figures
balance to 0.20 %, and the carbonate that charge balance *demands* sits within
4.1 mol/m³ of the curve actually read. Before the correction, Fig. 2's solid bulk
was short of anion charge by a factor of **5.25** — a factor of five sitting in
plain sight, invisible because the test that would reveal it was unbuildable.

So: **before trusting a multi-curve extraction, count the curves against a
conservation law the figure must satisfy.** If the law needs a species you did
not extract, that absence is the finding. And check the frame — a curve can exit
through the top of the box, and a pair can merge at the edge where you are most
likely to read them.

Worth noting what did *not* move: the α = 0.2 identification and the scalar `s`
were unchanged, because the two mis-read curves were *promoted* carbonates, which
enter neither. A correction that changes nothing downstream is still worth
making — it is what licenses the claim that nothing downstream changed.

### An old page's errors travel into the new page that reuses it

`F1.3` lifted a sentence from the published `F1.4` — the reuse `AGENTS.md`
encourages — and its verifier, reading Fig. 6 for itself, found the band wrong at
both ends (`U_df` 0.02–0.035 m/s against a measured 0.0161/0.0230/0.0269). It had
survived `F1.4`'s own review because it reads as background prose rather than as a
result, which is exactly where a wrong number hides.

So verifying a new page audits the pages it borrows from, for free — but only if
the verifier re-reads the *source*, not the sibling page. Say so in the dispatch.
And when a finding lands on a published page, fix it there in its own commit: it
is live, and the fix is cheap.

### Two independent readings of the same paper can validate each other

`J3.4` is the pattern. Its open-circuit expression (eq. 16) was read off a 600 dpi
render; the *dashed curve on Figure 2 is that same quantity*, and had already been
digitised and maintainer-reviewed. The two agree to **3.3 mV over 237 points**
against a digitisation good to ~3 mV. Neither reading informed the other, so each
validates the other — the page-image transcription and the curve trace are
independent witnesses.

Look for this whenever a figure plots something the text also states in closed
form. It converts two separately doubtful extractions into one confident result,
and costs nothing beyond noticing the overlap.

**When the paper does not print what the model needs, invert what it does print.**
`J3.4`'s Appendix A says the conductivity "was fit to a third-order polynomial"
and gives no coefficients. Rather than guess, the agent inverted the paper's own
eqs. 28 and 29, which both contain (1/kappa + 1/sigma) — two published numbers
yielding two estimates 19 % apart, stated openly as a reconstruction. That is the
line between reconstruction and fabrication: every input traceable to something
printed.

### `G1.8` — resolved: the legend is offset by one row against the curves

*This section replaces an earlier one titled "G1.8 is blocked". Every number in
it was superseded on 2026-08-02; the old text is gone rather than annotated,
because a stale finding sitting in the handoff is how the next agent re-derives
the wrong thing.*

Herskowitz & Smith's Figure 6 reproduced for exactly one of its four curves, and
the other three missed by **1.69×, 4.14× and 2.52×** at φ = 10. The earlier
conclusion — that nothing reconciled it, that the required α_gs showed no
pattern, and that the required *f*_e correction "grows with L_m" — was wrong on
the last point (it is **non-monotone**) and wrong in framing the whole thing as
irreducible.

**Ask which L_m each drawn curve *is*, instead of how wrong the model is.**
Solving the chain for the liquid rate that puts each printed curve where it sits
gives **0.508, 2.046, 7.061, 10.18** against the legend's 0.50, 1.0, 2.0, 7.0.
Three of the four are the paper's own printed values attached to the wrong rows.
Read that way the deviations fall from +69/+314/+152 % to **+2/+2/+6 %**, and an
independent re-digitisation from scratch put them at 0.14 %, 2.35 %, 0.59 % of
the printed values.

Two lessons generalise beyond this case.

**When a model misses a family of curves, invert it.** The residuals were not
noise and not a single bad constant — they were the right answers against the
wrong labels, which only shows up if you solve for the parameter instead of
scoring the fit.

**A separable model turns a figure's own spacing into evidence.** χ = g(L_m)·h(φ)
separates, so a gap between curves in decades is log₁₀ g(L_i)/g(L_j) whatever
h(φ) is. The printed reading forces d log g/d log L_m = −1.380, −2.024, −0.725 —
not constant, so no power law, and not even monotone. The shifted reading forces
−0.690, −1.120, −2.546 against the model's −0.684, −1.122, −2.414. This needs no
fitting and survives errors in the ordinate calibration.

The maintainer was asked with the figure inlined in a private decision artifact
and four options, and chose the shift. **That is the pattern to repeat** — see
[[gallery-decision-requests-need-the-artifact]]: name the PDF, show the figure,
state the alternatives and what each changes.

Two printed defects are recorded on the page: the legend offset, and page 8's
citation of "Table 2" for the spherical solution when Table 2 holds the
pressure-drop constants and the sphere row is in Table 1.

**And one caveat the page carries:** the residual log–log slope gap does *not*
go away under the shift, because the chain factorises into an L_m prefactor times
a φ shape — all four L_m give an identical slope, so reassignment moves it by
exactly zero. It is +0.007 to +0.046 depending on the abscissa calibration, and
that abscissa is the weak axis. Unexplained, and stated as such.

### Working the catalogue with parallel agents — the operating procedure

**Never block on the maintainer.** A case that needs input is *parked*, not
waited on: write its `resume:` block, and immediately dispatch the next case.
Keep about five builders live continuously. The maintainer reads the dashboards
when they choose; do not ping them per case and do not idle.

**The loop.**

1. Pick the next case (policy below) and dispatch a builder with
   `docs/agent-brief.md`.
2. When it reports `ready`, dispatch a **verifier** on the staged page before
   publishing. Findings go back to the builder or get fixed inline.
3. Integrate: move `queue_cases/<ID>/page/` to `pages/<ID>-<slug>/`, splice the
   agent's `models_entry.yaml` into `models.yaml`, flip `meta.yaml` to
   `published`, run `check_metadata.py`, `run_pages.py --changed`,
   `check_agreement.py`, commit, push.
4. When it reports `needs-input` or `needs-paper`, it is already parked — confirm
   the `resume:` block is complete, regenerate the dashboards, republish, move on.
5. Refill to five live.

**Dispatch policy, highest first.** Throughput is set by how rarely a case needs
the maintainer, so prefer cases that will not:

1. paper on disk **and** validation likely to be a table, appendix or stated
   result rather than a figure — these run end to end with no human involvement;
2. paper on disk, figure-validated;
3. paper reachable open access (`scripts/find_papers.py --fetch`);
4. anything else — these mostly become `needs-paper`, which is fine but yields
   no page.

Never dispatch two agents onto the same case, and never onto a case whose
`covered_by` names another.

**Draining parked cases.** When answers arrive, read each entry's `resume:`
block: `established` says what not to redo, `answer_changes` says exactly what
the answer alters, `files_to_touch` says where. A parked case should never need
its extraction or its validation repeated.

**Decisions that belong to no case go in `docs/standing-decisions.yaml`.** A
per-case blocker rides along in `queue_cases/<ID>.yaml` and reaches the dashboard
from there, but a repo-wide decision — a history rewrite, a permission grant —
has nowhere to live and used to survive only by being retyped into the next
session's prompt. Entries in that file render at the top of the decisions
dashboard with their own answer boxes. Delete one when it is settled and record
what was decided here.

**Integration gotchas, all seen at least once.** `slug` and `title` in
`meta.yaml` must match `models.yaml` exactly. An `id` may already exist as
`planned` — upgrade it in place rather than appending, or `check_metadata`
reports a duplicate. New dependencies must reach `requirements.txt` or CI breaks
on the next machine (`sympy` for `J4.8`).

**NEVER `git add -A` after an agent has run.** Review overlays are drawn on the
source page image, so they *are* the copyrighted figure; three reached this
public repo before being caught. `queue_cases/*/review/*.png` is git-ignored now,
but look at what you are staging.

**And check what *generated* files inline, not just what agents write.** The
git-ignore on the PNGs was not enough. `docs/dashboards/needs-input.html` is
tracked, and `scripts/dashboards.py` inlined those same overlays into it as
base64 — three in `7a2ed33`, then eight in each of `4d24b87`, `d0ae665`,
`1359629`, `cdcfd4b`. Ignoring a file and then embedding its bytes in a tracked
file publishes it just the same, and the deletion commit `e9a3cc8` did not touch
this route at all. Found 2026-07-31.

Fixed by splitting the build: `python scripts/dashboards.py` writes the tracked
pages with **no images** (each overlay is named instead), and `--with-images`
additionally writes `docs/dashboards/private/`, which is git-ignored and is what
gets published as the maintainer artifact. The tracked page went 903 kB → 23 kB.
The history still holds the old copies — that is a standing decision, below.
*The general rule: `git check-ignore` proves nothing about a generated artefact.
Grep the tracked output for `data:`, and check the size.*

**Three acquisition failure modes worth knowing.**

*A DOI resolved from a terse citation is usually wrong.* CrossRef returns
something confident for any query: "Carman (1937)" matched a 2025 paper citing
Kozeny–Carman. Year agreement is required, and auto-resolved DOIs are marked
unverified on the dashboard.

*Filenames carry no metadata.* Half the PDFs on disk are named by publisher PII
(`i260028a001.pdf`). `docs/papers-on-disk.yaml` maps catalogue ID to filename by
hand and is consulted first; add a line whenever a PDF arrives.

*Do not ask for papers a case does not need.* `A1.2`/`A1.3`/`A1.4` sat on the
papers dashboard while the catalogue asks for them as one comparison page — which
`A1.1` is. Check `covered_by` before reporting `needs-paper`.

*Open access will not carry section A.* Most 1937–75 classics predate DOIs
entirely. Expect that list to be irreducible and old.

### Working the catalogue with parallel agents

`queue_cases/` holds one YAML per catalogued case; `scripts/case_queue.py` and
`scripts/dashboards.py` drive it, and `scripts/find_papers.py` does acquisition.
Agents build into `queue_cases/<ID>/page/` and touch nothing shared, so several
run at once; an integrator merges finished pages into `pages/`, adds the
`models.yaml` entry, and runs the gates.

**Three things that went wrong in the first batch, all worth avoiding again.**

**Never `git add -A` after an agent runs.** Review overlays are drawn on the
source page image, so they *are* the copyrighted figure — three were committed to
this public repo before being caught, contradicting the redistribution basis every
sidecar states. `queue_cases/*/review/*.png` is now git-ignored, but the gate is
human: look at what an agent produced before staging it.

**A DOI resolved from a terse citation is usually wrong.** CrossRef returns
something confident for any query: "Carman (1937)" matched a 2025 paper citing
Kozeny–Carman, and "Danckwerts" matched a 1968 re-derivation. Title-word overlap
is not enough — require the publication year to agree, and mark anything
auto-resolved as unverified so a wrong DOI never sends the maintainer after the
wrong paper.

**Filenames carry no metadata.** Half the PDFs on disk are named by publisher PII
(`i260028a001.pdf`) or an export id, so automatic matching missed them and one
case was reported as needing a paper the maintainer had already supplied.
`docs/papers-on-disk.yaml` maps catalogue ID to filename by hand and is consulted
first; add a line whenever a PDF arrives.

**Open access will not carry section A.** Of the first cases checked, most 1937–75
classics have neither DOI nor open copy — they predate DOIs. Expect the papers
list to be dominated by old classics and to work well only for recent sections.

### Verify an equation read off a page image before building on it

`B3.1` is the template. Its two governing equations were read from 600 dpi
renders because the scan's text layer mangles them (`theta_B` as `0B`, Eq. 6's
exponents dropped). A page-image reading is a transcription and needs checking
like any other. Two checks settled it:

- **Endpoint identities that only hold for the right coefficients.** Eq. 6 must
  give exactly 0 at *r*/*R* = 1 and exactly 1 at *r*/*R* = 0 for *any* parameter
  values, which requires its numerator to collapse to its denominator term by
  term. A mis-read coefficient breaks this immediately.
- **An independent derivation.** Integrating the moving boundary from Eq. 5's
  three resistances, without looking at Eq. 6, reproduced it to 6.9e-16 over six
  decades — and recovered the factor 3 on the film term and the 12 in
  *k*_d1 = 12D/D_p on the way.

Look for both before writing the notebook. This is the same instinct as the
`C2.1` Table 5/6 round trip, generalised: *a transcription you can only read once
should be checked against something you can derive.*

**And refuse to invent what the paper omits.** `B3.1` cannot produce absolute
burnout times, because Parker & Hottel's correlation is printed for the specific
combustion rate `K_c` and the unit conversion to `k_c1` is not given. The page
works entirely in the dimensionless groups instead, which costs nothing, and says
so. Do not reconstruct a missing unit conversion by inference.

### The batched figure review, 2026-07-29 — read this before digitising anything

Five figures were digitised, put to the maintainer as **one** review artifact
(source figure ⇄ overlay toggle, closed questions plus a free-text box per
figure), and all five came back answered in a single pass. That is the workflow
to repeat: batching costs nothing extra and turns five review round-trips into
one. The artifact is private and no page image enters the repo.

What came back, and what it changed:

- **`H1.7` Wijmans & Baker Fig. 5** — all correct. Page now built and published.
- **`F2.3` Maretto & Krishna Fig. 2** — *"On the eps_s=0.35 line you detect a few
  circles and diamonds. These are just squares almost on top of other squares."*
  Fixed, see the lesson below. Also: report the unresolved cluster flagged rather
  than drop it, and the two circles floating above the top curve are real data.
- **`J3.4` Doyle–Fuller–Newman Fig. 2** — confirmed **no experimental data at
  all**; reference-solution page only. One print speck removed.
- **`G1.8`** — **switch to Figure 6, do not digitise Figure 2.** Fig. 2 needs the
  identity of ~45 overlapping markers and the correlation under test depends on
  particle diameter, so the F1.4 shortcut does not apply. Fig. 6 is the model.
- **`B3.1`** — agreed it is analytic; **no figure needed, no review gate.**

Staged extractions with their review verdicts are in `docs/staged-data/`.

**Series identity can come from position instead of shape.** This is the F1.4
lesson in a second form and it is the reusable one. On `F2.3` the template
matcher picked the wrong *shape* wherever markers overlapped — two overlapping
squares read as a circle. But the three series each follow their own curve and
never cross, so identity was reassigned by fitting `eps = a + b·log U + c·log²U`
per series and moving every marker to the nearest curve, iterated to a fixed
point. That moved exactly the 10 markers the reviewer had described, without
touching the detector. **Ask what carries the series identity in the figure —
shape, position, or a curve — before trying to improve shape recognition.**

**Check whether the figure can resolve what you are testing.** On `H1.7` the
rejection panel appears to validate the model; closing the model shows it
predicts a 0.31 percentage-point rise, which is 4 px on a figure whose curve is
6 px thick. Agreement there is not evidence. Compute the predicted effect in
*pixels* before quoting an agreement.

**Reusable extraction code** lives in the session scratchpad, not the repo:
`markers.py` (shape-template matching pursuit, open and filled), `curves.py`
(column tracing with an order-preserving DP so adjacent curves cannot swap), and
`overlay.py`. Worth promoting into `scripts/` when a fourth figure needs them.

### What `F1.4` settled, and what it did not

Figure 11 is digitised (63 markers) and the page is live. Three things are worth
carrying forward.

**The figure was salvaged by dropping labels, not by fixing the classifier.**
The four series differ only by marker shape, and shape recognition failed in the
dense band — the maintainer review confirmed it. Rather than curate a subset,
every row except the SF₆ group is `unassigned`. This cost nothing, because
**Eq. 19 contains no gas-density term**: the correlation can be tested against
positions alone, using all 63 points. Labels are needed only for the
gas-independence claim, and the SF₆ group alone carries that, being the density
extreme. *Generalise this before spending hours on a classifier: ask which
columns the model you are testing actually reads.*

**A group comparison can be confounded even when the groups look clean.** The
first version of the gas-independence test compared the SF₆ points' bias against
the rest — until a check showed the two groups occupy *disjoint* velocity
windows on this figure, every SF₆ point below 0.044 m/s and every other point
above 0.051. The difference measured velocity, not density. Replaced by an
extrapolation test: fit a free power law on the lighter gases only, predict the
SF₆ points it never saw. Before reporting a between-group difference, check the
groups overlap in every other variable.

**Deviation direction is not cosmetic.** The first draft computed the two
correlations' deviations in opposite senses — measured/model for one,
model/measured for the other. At 14 % scatter the reciprocal differs materially:
the mean moved 13.3 → 13.8 %, the bias −0.2 → +2.8 %, and the headline
gas-density number 3.5 → 5.0 % (before that test was replaced outright). Fix a
single convention,
(model − measured)/measured, state it on the page, and use it everywhere.

Figures 7 (column diameter) and 9 (liquid properties) are still undigitised, so
the diameter and liquid-property independence claims remain untested.

## Hard-won lessons — read before building

The first two are also stored as user memories and will load automatically.

- **Never fabricate data.** If a dataset cannot be obtained, mark it
  `status: placeholder`, keep the page `status: planned`, and say so on the page.
  Re-simulating the model and presenting it as data is circular and would
  discredit the gallery.
- **Reported numbers must be deterministic.** On `B1.1` the fold was first
  located on a warm-start continuation curve; CI re-executed on another machine
  and got η_ignited 44.45 against 36.15 locally. Locate features on a smooth
  deterministic reference (analytical, or a shooting solve) and rank turning
  points by prominence.
- **Validation catches what inspection does not.** Three defects in `B1.1` were
  found by the monotonicity and closure checks, not by reading the code: a sign
  error, a sweep that missed a solution branch entirely, and a comparison that
  reported 89% deviation where the truth was 1e-5.
- **`check_agreement.py` failing is a signal, not a nuisance.** Update the
  baseline only when the model genuinely changed, and say why in the commit.
- **Look for a check the *paper* pays for.** `C2.1` gained three that cost
  nothing: Tables 5 and 6 are related by the Arrhenius form, so recomputing one
  from the other tests the page-image reading *and* the split reference
  temperature at once; the van 't Hoff slope of an equilibrium correlation is a
  reaction enthalpy, so it can be checked against the paper's own table; and two
  figures plotting the same runs must pair up point-for-point across
  independently fitted axes. Look for these before writing the notebook — they
  are worth more than any amount of code review.
- **Two independent methods beat one method plus a tolerance.** On `D2.2` the
  critical inlet pressure is computed twice: by bisection on the pymrm reactor,
  and by an adaptive quadrature of the trajectory in the *p*–*T* phase plane
  that never forms the reactor grid. They agree to 0.18 % across the whole
  operating range, and that number is the only one on the page that does not
  involve the paper. Look for a second route to the same quantity — it is often
  cheap, and it catches discretisation error that grid refinement alone will
  flatter.
- **Marker extraction is per-figure, not a recipe.** Morphological opening
  worked on Duncan & Toor because its markers are solid and much larger than the
  curves. Xu & Froment's are ~20 px glyphs on ~10 px curves and needed a
  different method: local ink density minus the largest value explainable by a
  locally straight structure (grey-scale opening with a long *line* element,
  maximised over six orientations). Take the maximum over orientations, or every
  steep near-origin curve section is flagged as a marker. Then audit every
  candidate visually at 600 dpi — on `C2.1` roughly a quarter of the automatic
  candidates were curvature artefacts with no glyph at the crosshair.

## Working commands

```bash
cd ~/Code/pymrm_suite/pymrm-gallery
source ~/Code/pymrm_suite/.venv/bin/activate

python scripts/check_metadata.py            # metadata + provenance validation
python scripts/check_metadata.py --report-missing   # catalog IDs not yet in models.yaml
python scripts/run_pages.py                 # execute every page notebook (slow: ~minutes)
python scripts/check_agreement.py           # metric regression check

export QUARTO_PYTHON=$(which python)
quarto render                               # build the site into _site/
```

Pages are generated from builder scripts rather than hand-edited JSON.
`pages/C2.1-xu-froment-smr/build_page.py` is the reference example: a list of
`nbformat` markdown/code cells written to `index.ipynb`, so prose and code stay
in one reviewable file. Nine required sections, listed in
[`AGENTS.md`](../AGENTS.md).
