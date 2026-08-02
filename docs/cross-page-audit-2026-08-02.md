# Cross-page dataset audit, 2026-08-02

Sweep of the seven published pages (other than `A1.6`, already fixed) that call
`gallery_utils.load_data(..., page=<another page>)`, for the defect class where a
borrowing page treats a cross-page dataset as a bare table of numbers rather than
as a page that has already established findings about those numbers.

Every claim below was demonstrated by executing code or by reading the source
paper. Suspicions that could not be demonstrated are marked PLAUSIBLE and kept
separate.

## The borrow map

| borrowing page | source page | dataset(s) | columns used | what the borrower concludes from them |
|---|---|---|---|---|
| `A4.4` | `A4.2` | `krishna-wesselingh-1997-worked-examples.csv` | `value` for `ideal/D12,D13,D23` | the three H₂/N₂/CO₂ pair diffusivities, used only for the dimensional illustration (pore size at Kn = 1, the transition figure, the p- and d₀-scaling test) |
| `A4.3` | `A4.2` | same file, same three rows | same | the same three diffusivities, used throughout the ternary dusty-gas solves — the uphill-diffusion threshold, the viscous-relief numbers, the free-molecule and Maxwell–Stefan limits |
| `A4.2` | `A4.9` | `duncan-toor-1962-run1.csv` | `time_h, bulb, species, x` | the 28 digitised bulb compositions; the t = 0 rows set the initial condition and the rest are the page's only experimental comparison (0.59 mole % MS vs 3.14 % Wilke) |
| `F1.3` | `F1.4` | `krishna-ellenberger-1996-fig11.csv`, `-parameters.csv` | `U_minus_Udf, eps_b, gas`; column diameter, Reilly B, gas densities, Table 3 deviations | Wilkinson's eqs. 1–4 overpredict ε_b by +63.9 % (bias +46.8 to +79.4 % over the four gases); the error is a level error, not a slope error |
| `E1.2` | `E2.1` | `kunii_levenspiel_1968_appendix_values.csv` | 8 rows: `A/B ub_minus_u0`, `A/C K_bc`, `B H_bc`, `C u_br`, `C u_b`, `C gamma_c` | the Davidson potential-flow solve reproduces K&L's printed appendix values (worst 1.06 %) |
| `J3.1` | `J3.4` | `doyle-fuller-newman-1993-parameters.csv`, `-stated-results.csv` | full parameter set; `delta`, `nu` | the two electrodes of one cell sit ~3 decades apart in |η_s|; ν²/δ from Eq. 30 is −15.7 % against the printed groups |
| `J3.5` | `J3.4` | same two files | full parameter set; `S_c`, `S_s`, `u_at_sharp_drop` | the SPMe's asymptotic conditions are marginal at I = 10 and violated at I = 20; SPMe RMS error vs the DFN |

Group structure: `A4.2` is the source for `A4.3` and `A4.4` *and* a borrower from
`A4.9`; `J3.4` is the source for both `J3.1` and `J3.5`.

## Findings, ranked

### 1. `J3.1` — the I = 20 collapse utilisation contradicts the dataset it loads. CONFIRMED

`J3.1` states, in *The result: the two electrodes of the same cell are three
decades apart*:

> At $I=20$ the cell is transport-limited and collapses at $u=0.40$; the cathode
> peak overpotential rises to 3.7 mV there … That is the largest kinetic
> overpotential anywhere in any of these runs.

The page loads `doyle-fuller-newman-1993-stated-results.csv`, which carries
`u_at_sharp_drop, 0.30, text journal page 1529, I = 20 A/m2 where the cell
potential drops sharply`, and `doyle-fuller-newman-1993-parameters.csv`, which
carries `V_cutoff, 1.7 V`. It uses neither. `Cell.march` hard-codes
`v_stop=1.55`, 150 mV below the paper's own cutoff, and the page reports the
march endpoint as a result.

Verified in the paper (Doyle, Fuller & Newman 1993, p. 1529, read from the PDF
text): *"at a rate of 20 A/m² the cell potential drops sharply when about 30% of
the cathode material is utilized … A typical cutoff voltage is about 1.7 V;
beyond this value the cell is severely polarized."* `J3.4` independently
reproduces this: its published Figure 2 curve for I = 20 is drawn only to
u = 0.299 and its own 1.9 V crossing is at u = 0.264.

Measured by re-executing `J3.1`'s own cells 1–12 and marching at I = 20:

| quantity | value |
|---|---|
| u at V = 1.9 V | **0.2638** (= `J3.4`'s 0.264) |
| u at V = 1.7 V, the paper's cutoff | **0.3707** |
| u at the march endpoint, V = 1.473 V | 0.4033 — the page's "0.40" |
| max \|η_s\| up to the 1.7 V cutoff | **1.679 mV** |
| max \|η_s\| over the full run | 3.705 mV — the page's "3.7 mV" |

At I = 10 the same run gives u = 0.8312 at the 1.7 V cutoff — exactly `J3.4`'s
0.831 and the paper's stated 0.84 — so the model is right and only the readout
convention is wrong. `J3.4` marches with the same `v_stop=1.55` but then reads
its results at `P["V_cutoff"]` via `crossing(...)`; `J3.1` inherited the march
and not the readout discipline.

**Failure scenario.** A reader takes "3.7 mV at u = 0.40" as the worst kinetic
overpotential inside the published operating envelope — the section is titled
*The published cell as an operating envelope*. Inside the envelope the worst case
is 1.68 mV, a factor 2.2 smaller, and the collapse is at u ≈ 0.30–0.37, not 0.40.
The page's qualitative conclusion (both electrodes deep inside the linear window)
survives; the two numbers do not.

Nothing else on the page is affected: the headline overpotentials at I = 5 and
I = 10 are taken under a `u < 0.80` mask and are unchanged, and
`agreement.json` records only the I = 10 values.

### 2. `J3.5` — Reuse advertises `bv_invert` without the condition that makes it exact. CONFIRMED

`J3.5`'s Reuse section:

> **The closed-form Butler–Volmer inversion.** `bv_invert` solves Doyle's Eq. 17
> for the overpotential in closed form (a quadratic in $e^{\eta'}$ …). Any model
> with that kinetics — including the J3.4 DFN itself, in a preconditioner — can
> reuse it.

`bv_invert` forms `z = (B + sqrt(B² + 4 c_s (c_T − c_s)))/(2 c_s)` and returns
`U' + (RT/αF) ln z`. That is a quadratic in `exp(η')` only when the two
exponentials are `exp(+αη')` and `exp(−αη')` with the same α — i.e. when
α_a = α_c, equivalently α_a + α_c = 1 with equal transfer coefficients. `J3.5`
never states the restriction, and its own α_a = α_c = 0.5 hides it.

`J3.1` — the sibling page built on the same `J3.4` dataset — has already
identified this and says so in *its* Reuse section: *"note that it is exact only
under the α_a + α_c = 1 condition derived here, which that page does not
state."* The finding exists in the repository and has not travelled back to the
page it is about.

**Failure scenario.** A reader with asymmetric transfer coefficients (α_a = 0.3,
α_c = 0.7 — a routine case, and the one `J3.1`'s Check 2 exercises) copies
`bv_invert` on `J3.5`'s explicit invitation and gets a silently wrong
overpotential: the quadratic root is not the solution of Eq. 17 for that α pair,
and there is no error, only a wrong number.

### 3. `A4.2` — the notebook still calls the pair diffusivities "measured", contradicting its own corrected sidecar and the paper. CONFIRMED

Two places in `index.ipynb`:

- markdown: *"D₁₂ = 8.33, D₁₃ = 6.8, D₂₃ = 1.68 × 10⁻⁵ m² s⁻¹ ('from the kinetic
  gas theory' — numerically identical to the values Duncan & Toor **measured** at
  35.2 °C and 1 atm)"*;
- code comment: *"# Identical to Duncan & Toor's **measured** 0.833 / 0.680 /
  0.168 cm²/s."*

The page's own sidecar, `krishna-wesselingh-1997-worked-examples.meta.yaml`, was
corrected on 2026-08-01 and now reads: *"Corrected 2026-08-01: an earlier version
of this file said Duncan & Toor 'measured' them, which neither paper
establishes."* The correction landed on the sidecar and not on the notebook.

The paper settles it. Duncan & Toor (1962), p. 40, verbatim: *"The binary
diffusion coefficients used were the best experimental values available (16)
corrected to the thermostat temperature by the method suggested by Hirschfelder,
Curtis, and Bird (7). The values of the diffusion coefficients at 35.2 °C and
1 atm were nitrogen-hydrogen, 0.833 sq. cm./sec.; …"* They are literature values
temperature-corrected by an HCB procedure — neither measured by Duncan & Toor
nor, strictly, "from the kinetic gas theory" as Krishna & Wesselingh describe
them. `A4.9`, the page that owns the dataset, never claims they were measured.

**Failure scenario.** A reader building a third page treats 0.833/0.680/0.168 as
a tier-4 measurement of the H₂/N₂/CO₂ pairs at 35.2 °C, cites Duncan & Toor for
it, and states an experimental provenance the literature does not support. The
whole reason `A4.4` and `A4.3` phrase this carefully is that the sidecar told
them to; the notebook says the opposite.

### 4. `A4.2` — 0.45 mole % and 2.6 mole % attributed to the wrong paper. CONFIRMED

`A4.2` prints, beside its experimental comparison:

> For scale: **the review's own predictions** deviate 0.45 mole %, the
> experimental error is 2.6 mole % …

On `A4.2`, "the review" means Krishna & Wesselingh (1997) throughout — the phrase
is used ~20 times and always for that paper. Both numbers are Duncan & Toor's:
*"For all ternary diffusion runs the average deviation of the predicted mole
fractions by the Maxwell-Stefan equations from the measured mole fractions is
0.45 mole %. This is well within the expected experimental error of 2.6 mole %."*
(Duncan & Toor 1962, p. 40.) Neither appears in the K&W review. `A4.9`, the
source page, attributes them correctly: *"The paper reports 0.45 mole % average
deviation of its own Maxwell-Stefan predictions … and quotes an experimental
error of 2.6 mole %."*

**Failure scenario.** A reader goes to Krishna & Wesselingh for the 0.45 mole %
benchmark, does not find it, and either drops the comparison or re-derives it.
The number itself is right; only its owner is wrong.

### 5. `A4.4` — "no validated number depends on the temperature or pressure" is not true as stated. CONFIRMED (overstatement, conclusion survives)

`A4.4` says of the 35.2 °C and 1 atm it inherits through `A4.9`'s sidecar: *"they
enter only through a mean molecular speed and a total concentration — **no
validated number on this page depends on either**"*, and repeats it in the code
comment.

Break test, re-executing all of `A4.4` with the value substituted and comparing
the 40 metrics in `agreement.json`:

| injection | metrics that move | largest move |
|---|---|---|
| `T_REV` 308.35 → 400 K | 2 of 40 | `scaling_slope_dev_knudsen` 0.00730 → 0.00830 (+14 %); `scaling_slope_dev_bulk` 0.0626 → 0.0558 (−11 %) |
| borrowed D₁₂/D₁₃/D₂₃ × 1.3 | 2 of 40 | `scaling_slope_dev_knudsen` 0.00730 → 0.00564 (−23 %); `scaling_slope_dev_bulk` 0.0626 → 0.0787 (+26 %) |

The two that move are the page's check against the four exponents the review
states in words on p. 887, measured as log–log slopes over a sweep whose window
is set by the dimensional diffusivities. Both stay far below 0.1, so the
conclusion — the review's stated exponents are reproduced — is unaffected. The
other 38 metrics are bit-identical, which is a genuinely strong result and worth
publishing as such. The sentence should say "no *conclusion* depends on either,
and two exponent-recovery deviations shift by ~10–25 % without changing the
verdict", not "no validated number".

### 6. `F1.3` — the two Figure 6 panels are labelled opposite ways in prose and in code. CONFIRMED

*The data*, item 4: *"three (ρ_G, ε_df) pairs from **panel (b)** and three
(ρ_G, U_df) pairs from **panel (a)**"*.

The code cell immediately below: *"The caption reads 'Influence of gas density on
the (a) dense-phase gas voidage, and (b) superficial gas velocity through the
dense phase', so eps_df is **panel (a)** and U_df is **panel (b)**."*

The caption is quoted correctly (verified against the paper). The printed figure
carries no (a)/(b) marks on the panels themselves, and the velocity panel is
drawn *above* the voidage panel, which is presumably where the prose reading came
from. Either reading is defensible from an ambiguous figure; carrying both on one
page is not.

The six numbers themselves are correct. Read off a 300 dpi render of journal page
2631: the lower panel gives ε_df ≈ 0.09 / 0.14 / 0.19 at ρ_G = 0.18 / 1.30 / 1.83
and the upper gives U_df ≈ 0.016 / 0.023 / 0.027 m/s, matching the page's
0.089/0.139/0.192 and 0.0161/0.0230/0.0269. `F1.4` carries the same six values and
records that an earlier version of *itself* quoted the wrong bands — the
inheritance defect that started this class is properly closed on both pages.

### 7. `E1.2` — a printed slip recorded on `E2.1` is reproduced without the note. CONFIRMED (immaterial)

`E1.2`'s provenance table prints, as one of the rows it says it re-read on a
600 dpi render: *"`C u_br`, `C u_b` | u_br = 0.711(980×3.7)^{1/2} = 42.8;
u_b = 13.2 − 2.1 + 42.8 = 53.9"*, and uses `u0 = 13.2` in the recomputation.

`E2.1`'s sidecar records that appendix C prints *"u₀ = (6.6 + 9.9 + 13.2 + 20)/5
= 13.2 cm/s — four terms divided by five"*, notes that 13.2 is nonetheless the
value the authors used throughout and is confirmed by u_b = 53.9, and refuses to
infer the missing term. Verified in the paper (Kunii & Levenspiel 1968, appendix
C): the expression is printed exactly as `(6.6 + 9.9 + 13.2 + 20)/5 = 13.2
cm./sec.`

`E1.2` carries none of that. Nothing on the page depends on u₀ beyond the one
reproduction row, so the consequence is small — but the page asserts it re-read
those rows and the slip is in the row it re-read.

`E1.2` handles the *other* `E2.1` slip well: it explicitly tests whether appendix
B's verdict survives the ε_mf = 0.447 implied by the printed 8.70, and reports
that it does. That is the behaviour this audit is looking for.

### 8. Provenance-tier convention differs between two borrowing pages. CONFIRMED (catalogue metadata)

`models.yaml` gives one `data.tier` per page, and the two pages whose only
experimental data is borrowed treat it differently:

- `F1.3` — own dataset is a table, borrowed dataset is `F1.4`'s digitised Figure
  11 → `data: {tier: 4, method: digitised}`;
- `A4.2` — own dataset is the K&W worked examples (tier 6), borrowed dataset is
  `A4.9`'s digitised Figure 2 (tier 4, and the page's headline 0.59 mole %
  comparison) → `data: {tier: 6, method: table}`.

Both pages' prose is honest — `A4.2`'s `meta.yaml` caveats say *"The page's own
dataset is tier 6 … The experimental comparison uses the tier-4 digitised Duncan
& Toor dataset published with page A4.9, loaded cross-page"*. Only the
machine-readable field disagrees. A reader filtering `models.yaml` for
`data.tier: 4` to find pages validated against measurement gets `F1.3` and misses
`A4.2`.

### 9. `J3.4` and `J3.5` retype δ = 1.95 instead of reading it from the CSV they load. CONFIRMED (minor)

Both define `def kappa_from_delta(delta=1.95, I=10.0)` and call it with the
default, although `doyle-fuller-newman-1993-stated-results.csv` — loaded on both
pages — carries `delta, 1.95`. κ is the reconstructed, load-bearing input on
both. `J3.1` does it correctly (`lhs = S["nu"]**2 / S["delta"]`). If the
transcription of δ were ever corrected in the CSV, `J3.1` would follow and the
other two would not.

## Defect injections run

Every borrowed dataset was perturbed and the page re-executed end to end, with
`report_agreement` stubbed so no file in the repository was written.

| page | injection | result |
|---|---|---|
| `A4.4` | borrowed D × 1.3 | 2 of 40 metrics move (see finding 5); 38 bit-identical — the page's claim that the borrowed values are decorative is essentially true |
| `A4.4` | `T_REV` → 400 K | 2 of 40 move, 11–14 % |
| `A4.3` | borrowed D × 1.3 | **19 of 39** metrics move, including `uphill_threshold_nm_at_x_0p49` 212 → 276 nm and `viscous_relief_1um_pore` +26 % — the results genuinely rest on the borrowed numbers, as the page says |
| `E1.2` | all borrowed printed values × 1.1 | the 3 metrics that compare against printed values move (0.43 → 9.48 %, 1.06 → 10.05 %, 0.46 → 9.51 %); the 10 internal-identity metrics do not, correctly |
| `J3.1` | re-march at I = 20, read at the paper's cutoff | see finding 1 |

`F1.3`'s +63.9 % is a direct point-by-point deviation against the 63 borrowed
markers and moves with any perturbation of them by construction; no injection was
needed. `J3.5` publishes its own injection — it measures that the shared-κ
blindness is worth 0.8 mV and 9 % against the 9.8× a 10 % error in the
reduction's own ohmic coefficient produces — which is the right way to do it.

## What each page needs

| page | verdict |
|---|---|
| `A4.3` | **Clean.** Reflects every `A4.2` finding that touches the borrowed rows (kinetic-theory phrasing, no printed T or p, the A4.9 provenance of 35.2 °C / 1 atm, tier 6), the borrowed values are load-bearing and stated to be, and the Reuse section contradicts nothing on the source page. No action. |
| `E1.2` | **Clean but for finding 7.** It reflects `E2.1`'s appendix-B slip and tests that its verdict survives it; add the one-line note about appendix C's u₀ arithmetic. |
| `F1.3` | **Clean but for finding 6.** It inherits `F1.4`'s gas-label caveat explicitly, obeys it, uses the same deviation convention, cross-checks its eq. 19 number against `F1.4`'s and labels the match a transcription cross-check rather than corroboration. Fix the panel (a)/(b) sentence. |
| `A4.4` | **Needs one fix** — finding 5, a wording change plus (ideally) publishing the two-metric sensitivity table, which is stronger evidence than the current absolute claim. |
| `J3.5` | **Needs one fix** — finding 2, the `bv_invert` restriction in Reuse and in the docstring. Finding 9 optionally. |
| `A4.2` | **Needs two fixes** — findings 3 and 4, both provenance sentences, both already settled by documents on disk. Not a send-back: every number on the page is right. |
| `J3.1` | **Needs a specific fix, and it touches a result** — finding 1. The sentence, the two numbers in it, and preferably the readout convention (report at `P["V_cutoff"]`, as `J3.4` does, rather than at the march endpoint). Not a send-back: the model agrees with `J3.4` to three digits once read at the right voltage. |

No page needs to be sent back.

## The one structural change

Add to `AGENTS.md`, under **Data rules — the ones that matter most**:

> **Loading another page's dataset means reading that page.** `load_data(...,
> page=<other>)` is not a file path — it is a claim that the other page's
> findings about those numbers carry over. Before using a cross-page dataset:
>
> 1. Read the source page's `data/*.meta.yaml`, its *The data* section and its
>    *Validation* and *Reuse* sections. List, in your own *The data* section,
>    every finding the source has established about the rows you use — flagged
>    rows, printed slips stored as printed, reconstructed inputs, review
>    verdicts, quantities the source derives from the same columns — and say for
>    each whether your page's use is affected. "Not affected" is a fine answer;
>    silence is not.
> 2. **If a number you are about to state also exists in a dataset you loaded,
>    print it beside yours and reconcile them.** Never retype a value that is a
>    row in a CSV you already read; read it from the CSV.
>
> This is cheap and it is the only reliable guard. Every cross-page defect found
> so far — `A1.6`'s inverted voidage against `A1.7`'s computed one, `F1.4`'s
> wrong band travelling into `F1.3`, `J3.1`'s u = 0.40 against the `u_at_sharp_
> drop = 0.30` sitting in the file it loaded — is rule 2 not being applied.

Rule 2 alone would have caught findings 1, 3, 4 and 9 and the original `A1.6`
defect; rule 1 catches 2 and 7. The verifier brief's item 4 ("Borrowed claims …
verify against the *source*, not the sibling page") already covers the reviewing
side; nothing currently instructs the *builder*, and the builder is the one
holding the dataset.
