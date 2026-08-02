# F1.3 — adversarial verification

Verifier pass, 2026-07-31. Everything below was established from the 600 dpi
renders of `~/papers/pymrm-gallery/AIChE_Journal-1996-Krishna.pdf` and from
independent re-computation. `docs/pdf-findings.md`, the F1.4 code and the
builder's account were not used as evidence for any equation or constant.

---

## VERDICT ON THE CENTRAL CLAIM

**The Eq. 2 collapse is a REAL property of the correlation as printed, and it is
Krishna & Ellenberger's own published finding — not a units error, not a
transcription error, not a quantity mismatch, and not an out-of-range evaluation
in any sense the source states.**

Five independent eliminations:

1. **Transcription.** Read myself on 600 dpi renders. Eq. 2 on journal p. 2627
   reads exactly
   `U_trans = eps_trans V_small ;  eps_trans = 0.5 exp(-193 rho_G^-0.61 mu_L^0.5 sigma^0.11)`
   — minus sign inside the exponential, 193, prefactor 0.5, and all three
   exponents as the page prints them. Eqs. 1, 3, 4 also verified character by
   character (including `(sigma^3 rho_L/(g mu_L^4))^-0.273` and `^-0.077`, the
   `(rho_L/rho_G)^0.03` and `^0.077`, and the `0.757`). Reilly Eq. 8 on p. 2629
   verified including `B^1.5`, `0.59`, `1/2.84`, `rho_G^-0.04`, and
   **`sigma^0.12` sitting inside the radical** — the page's reading is right.
   Eq. 19 on p. 2632 is `0.268 D_T^-0.18 (U-U_df)^-0.22 (U-U_df)^{4/5}`, i.e.
   exponent 0.58 as used. Table 2 and Table 3 verified.

2. **Units — settled outright, and by stronger evidence than the page uses.**
   The paper's own **Notation list on p. 2634** prints
   `sigma = surface tension of liquid phase, N m^-1`,
   `mu_L = liquid viscosity, Pa s`,
   `rho_G = density of gaseous phase, kg m^-3`.
   The correlation is SI-bound by the source's own declaration. The page never
   cites this and argues the point from Figure 4 instead; see finding 9.

3. **The authors plot the same number.** Figure 4 draws K&E's own evaluation of
   Eq. 2. Tick-calibrated pixel reads:
   - panel air–tetradecane: Wilkinson line at 0.0024–0.0025 (page computes
     0.0028); Reilly line at 0.150 (page 0.147).
   - panel air–water: Wilkinson line sits on the 0.01 decade tick (page 0.0102).
   - Figure 6(b), a *different* figure, puts the Wilkinson eps_trans line at
     ~0.0032 at rho_G = 1.3 (page 0.0028) and its Fig. 6(a) U_trans line at
     ~0.0009 (page 0.00069).
   Three independent readings of the authors' own computed lines agree with the
   page's SI evaluation. A cgs or mPa·s reading would give 3e-153 (the page
   prints this itself).

4. **Quantity.** K&E write "For estimation purposes of eps_df and U_df we
   recommend the use of Eq. 8 assuming eps_df = eps_trans and U_df = U_trans"
   and "the Wilkinson correlation severely underpredicts the values of the
   voidage and gas velocity through the dense phase". They compare Eq. 2 against
   measured eps_df in Figs. 4 and 6 themselves. The page compares the same
   quantities the authors do. (The eps_small = eps_df(1-eps_b) vs eps_trans
   distinction is worth at most a factor 1.33 and is swamped.)

5. **Validity range.** No validity range for Eq. 2 appears anywhere in K&E.
   The page states this. K&E *do* state Reilly's ceiling — "the Reilly
   correlation has been developed for a data set with a maximum value of
   eps_trans = 0.32" — which the page quotes correctly.

The helium number does **not** inherit F1.4's label problem: measured eps_df at
helium is a labelled point on Figure 6, a different figure, read here as 0.089
(page's band 0.10-0.20). 0.089/1.5e-8 is a real comparison at helium's real
density.

---

## FINDINGS

### 1. HIGH — "Eq. 2 becomes sensible above ~14 kg/m3" and "the two closures are complementary" are not supported. CONFIRMED.

`rho_G_where_eq2_reaches_0p15 = 14.2` is the density at which Eq. 2 crosses a
**fixed** 0.15. But Figure 6(b) — the paper's own measurement — shows eps_df
*rising* with gas density: tick-calibrated reads give **0.089 / 0.139 / 0.192**
at rho_G = 0.18 / 1.3 / 1.83, i.e. eps_df ~ rho_G^0.30. Extrapolated to
14.2 kg/m3 the measured level is **0.32**, so Eq. 2 there would still be a factor
2.1 low. The crossing density is an artefact of holding the target constant.

Worse for the "complementary" claim: **Eq. 2 never exceeds Eq. 8 at any density
up to 1000 kg/m3.**

| rho_G | Eq. 2 | Eq. 8 | ratio |
|---|---|---|---|
| 0.18 | 1.5e-8 | 0.057 | 3.8e6 |
| 6.7 | 0.0743 | 0.324 | 4.4 |
| 14.2 | 0.150 | 0.464 | 3.1 |
| 25 | 0.213 | 0.609 | 2.9 |
| 70 | 0.317 | 0.998 | 3.1 |
| 1000 | 0.457 | 3.58 | 7.8 |

The ratio bottoms out near 2.9 and rises again; the curves never cross, which is
visible on the page's own left panel. Nothing on the page or in the source shows
Wilkinson becoming *correct* at any density — only that it stops being absurd.

Compounding this, `axhspan(0.10, 0.20)` is drawn across the **entire** 0.1–40
kg/m3 axis (and `axhspan(0.020, 0.035)` likewise on the right panel). Labelled
"measured dense-phase voidage", a full-width flat band visually asserts
density-independence that Figure 6 explicitly contradicts and that the paper's
title sentence ("increasing gas density significantly increases the dense-phase
gas voidage") denies.

*Failure scenario.* A reader takes the Reuse advice literally — "for a
pressurised column use Wilkinson's Eq. 2 above ~14 kg/m3" — and gets a
small-bubble holdup a factor of two low, with the page's own figure appearing to
endorse it.

*Fix.* (a) Clip the grey bands to the density range they were read over
(0.18–1.83, or at most to 6.7). (b) State that the measured eps_df rises with
rho_G, so 14 kg/m3 is a **lower bound** on where Eq. 2 could become adequate.
(c) Replace "complementary rather than rival, each valid where the other is not"
with what is actually demonstrated: Eq. 2 is below Eq. 8 at every density
examined, and the operative asymmetry is that Eq. 8 leaves its authors' stated
0.32 range at 6.5 kg/m3 and leaves physics at 70, while Eq. 2 stays bounded by
its 0.5 ceiling. This affects `meta.yaml adds:`, `models_entry.yaml
description:`, the notebook's "This is the central result" cell, the "What pymrm
adds" second paragraph, the Reuse section and the README headline.

### 2. HIGH — the C1 check "recovers eq. 4's printed exponent 0.757 as 0.7566" is circular, and is claimed as a transcription check. CONFIRMED by experiment.

`wilkinson_vb` contains `** 0.757`; the fitted decay slope simply reads it back.
Substituting other values into that one line:

| exponent typed into `wilkinson_vb` | slope "recovered" |
|---|---|
| 0.757 | 0.7566 |
| 0.657 | 0.6557 |
| 0.900 | 0.8999 |
| 0.400 | 0.3725 |

The check cannot detect a mis-read exponent. The residual 0.05 % is finite-excess
truncation (`rel = x/(1+x)` with `x ~ d^0.757`), not an agreement.

What the check *does* establish is real and worth keeping: that Eq. 4's first
term is exactly Eq. 3's V_small, and that the correlation is C1 across its own
regime boundary. But the page relies on it as one of "two internal identities
[that] carry the transcription" and as the stated defence against inheriting a
typo from K&E's reprint. That defence is half illusory — only the dimensional
identity survives, and it constrains Eqs. 3 and 4 only.

*Failure scenario.* K&E mis-set 0.757; the page reports "0.757 read correctly"
and the caveat section says the identities pass.

*Fix.* Reword in `meta.yaml validation:`, `meta.yaml agreement:`,
`meta.yaml caveats:`, `README.md` and notebook Validation §2 to: "confirms that
Eq. 4's first term is exactly Eq. 3's V_small and that the correlation is C1";
drop "which checks the transcription independently" and "0.757 is read
correctly". Rename the metric `continuity_exponent_recovered` or describe it as
a self-consistency value. (The genuine external transcription check is the
Figure 4/6 reproduction, plus — see finding 9 — the Notation list.)

### 3. MEDIUM — "the whole +64 % belongs to eq. 4" understates Eq. 3's share. CONFIRMED.

Over the 63 points, `V_small` (Eq. 3) is **46 % of V_b on average**, ranging
23 %–77 %, and it is 77 % at the lowest excess velocity. "V_b is 39 % low" is
fine — V_b is Eq. 4's subject — but "the whole +64 % belongs to eq. 4", "the
other being eq. 4's rise velocity, 39 % low" (`scope_decision`) and "the +64 % on
eps_b is eq. 4's" (`meta.yaml adds:`) exclude Eq. 3, which supplies about half of
the quantity and dominates at the low-velocity end. The "two independent defects,
separated" framing is really Eq. 2 vs Eqs. 3+4; the split was made between the
transition equation and everything else, and the third equation was folded
silently into the Eq. 4 bucket.

*Fix.* "Eqs. 3 and 4 jointly", or "V_b — whose first term is Eq. 3's V_small — is
39 % low".

### 4. MEDIUM — the headline +63.9 % carries an unstated +-16 pp from evaluating all 63 points at air density. CONFIRMED.

| all points evaluated as | bias | implied V_b shortfall |
|---|---|---|
| helium (0.18) | +46.8 % | 32 % |
| air (1.3) | **+63.9 %** | **39 %** |
| argon (1.83) | +67.1 % | 40 % |
| SF6 (6.7) | +79.4 % | 44 % |

Refusing to use the unreliable labels is the right call and matches F1.4. But
F1.4 drew no quantitative conclusion from the number, whereas F1.3 builds a
two-defect decomposition and a "+67 to +88 % end to end" on the 39 %. The
single-density assumption is stated but its size is not.

*Fix.* Print the four-gas span alongside the headline and note that the
attribution is +-16 pp on the bias / 32–44 % on V_b.

*Cleared, related:* the disjoint-window confound from `handoff.md` **is** present
in this dataset (the 12 SF6-labelled points span x = 0.0143–0.0436, all 51 others
0.0508–0.3662, completely disjoint), but F1.3 makes no group comparison —
`is_sf6` is computed and used only in a count. No test on the page depends on
marker shape. The scope claim on that point holds.

### 5. MEDIUM-LOW — "one to two decades below the measurement" contradicts the page's own helium headline. CONFIRMED.

Appears in the "central result" cell and again in "What pymrm adds" ("one to two
decades too low at their conditions"). True for air (1.7 decades) and argon
(1.4), but at helium — one of "their conditions", and the page's own headline —
it is 6.8 decades. Say "one to two decades at air and above, and seven at
helium".

### 6. LOW-MEDIUM — the "68 % cancels" is derived, not fitted, but is a near-tautology and rests on a model reference. PLAUSIBLE overstatement.

`canc = 1 - |a+b|/(|a|+|b|)` for opposite-signed a, b equals `2 min/(|a|+|b|)` —
here 2(0.0662)/0.1941 = 0.682, i.e. a restatement of the 0.128 : 0.066 ratio.
Nothing is fitted to Table 3 (Table 3's values are not inputs), so the
circularity charge fails. Two caveats do apply:

- the small-bubble error term uses **Reilly's** eps_df, a model, as the stand-in
  for the measured small-bubble holdup; the notebook says so, `meta.yaml adds:`
  does not;
- the notebook line "eq. 2 underpredicts eps_small by about as much as eq. 1
  overpredicts eps_b" is wrong by a factor of ~2 (0.128 vs 0.066).

The computed 3.6x vs printed 5.2x ratio is honestly labelled a pattern, not a
reproduction. Table 3 is correctly *not* reproduced anywhere.

### 7. LOW — the "measured U_df 0.020-0.035 m/s" band does not match Figure 6. CONFIRMED.

Tick-calibrated reads of Fig. 6(a): **0.0161 (He), 0.0230 (air), 0.0269 (Ar)**.
The stated band excludes the helium measurement and its upper edge is 30 % above
the largest measurement; a one-significant-figure read would be 0.02–0.03. The
eps_df band 0.10–0.20 against measured 0.089/0.139/0.192 is fair.

This band is inherited verbatim from the **published F1.4 page**
(`pages/F1.4-.../build_page.py`: "the paper's Fig. 6 measures eps_df ~ 0.10-0.20
and U_df ~ 0.02-0.035"), so the maintainer may want to correct it there too.
Nothing quantitative rests on it on either page.

### 8. LOW — "recomputed here independently and agrees exactly" overstates.

The A4.2 failure mode is **absent**: F1.3 imports nothing from F1.4 and defines
its own functions; and F1.4 genuinely never evaluates Eq. 2 (no `193`, no
exponential transition anywhere in its build script) — the scope decision is
factually correct. But F1.3 retypes the same Eqs. 3 and 4 from the same source
and applies them to the same 63 rows of the same CSV, so agreement to all printed
digits is arithmetically guaranteed and tests only that the two transcriptions
match. Call it a transcription cross-check rather than an independent
computation.

### 9. LOW / informational — a stronger, free units check was left on the table.

The Notation list on p. 2634 prints the units of sigma, mu_L and rho_G outright
(N/m, Pa·s, kg/m3). Citing it would settle the SI question by quotation rather
than by inference from a figure, and would strengthen Validation §1, which
currently asserts "Eq. 2 and Reilly's eq. 8 are both SI-bound" without a source.

### 10. LOW / informational — Figure 6 offers a third reproduction, and one mismatch worth knowing.

Fig. 6 draws K&E's own Wilkinson *and* Reilly transition curves against gas
density for the exact system the page uses. The page's Wilkinson values
reproduce (0.0032 vs 0.0028 for eps_trans, ~0.0009 vs 0.00069 for U_trans at
1.3 kg/m3), and the page's Reilly eps_trans reproduces to 1–6 % (0.093/0.130/
0.158 computed vs 0.082/0.125/0.158 read at rho_G = 0.5/1.0/1.5).

But the page's **Reilly U_trans is 5–29 % below K&E's own drawn Reilly curve**:
computed 0.0131 / 0.0198 / 0.0258 / 0.0299 at rho_G = 0.18 / 0.5 / 1.0 / 1.5
against 0.0159 / 0.0209 / 0.0308 / 0.0385 read from the figure, the gap growing
with density. The page's implementation is faithful to Eq. 8 exactly as printed
(I verified the `(1 - eps_trans)` factor on the render), so this is K&E's own
inconsistency, not the page's. It does mean the parenthetical "(which the paper
recommends and Fig. 6 supports)" is generous — Reilly's U_trans is also +-20 %
against the measured U_df points (0.0131 vs 0.0161 at He, 0.0284 vs 0.0230 at
air, 0.0320 vs 0.0269 at Ar). The abscissa-convention conclusion ("at most 15 %")
is robust to a 20 % shift in U_df, so nothing on the page turns on it.

---

## WHAT PASSED

- **Determinism.** Re-executed from a clean kernel in a scratch copy. Every
  printed number is byte-identical to the staged notebook, and `agreement.json`
  matches the re-run. No warm-start chains, no RNG.
- **No stale hard-coded numbers.** Every figure in the prose (0.0028, 1.5e-8,
  4/53/3.8 million, 4 nm/s, 14 kg/m3, 11 bar, 6.54, 70, 15, +64 %, 39 %,
  +67–88 %, 15 %, 68 %, 3.6x, 8 %, 0.621/0.636/0.58, 12.2 m/s, 10^-153, 0.7566)
  traced to the executed output.
- **Convention.** (model - measured)/measured used everywhere, stated on the
  page, matches F1.4.
- **Table 3 not reproduced**, and said so in four places. Correct: the 1,735 runs
  are unpublished.
- **Reilly's mirror-image behaviour (h)** verified analytically against the
  printed Eq. 8 with B = 3.85: eps_trans ∝ rho_G^0.48 gives 0.32 at 6.58 (page
  6.54), 1.0 at 70.5 (page 70), and the U_trans maximum at eps = 0.5 near 16.6
  shifted down by rho_G^-0.04 (page 15). All correct.
- **Structure.** Nine required sections in order; Colab cell first; `load_data`
  with explicit `page=`; `report_agreement` with 13 metrics; no Quarto syntax;
  no page image; sidecar present and accurate; no `- id: F1.3` in `models.yaml`
  so the append instruction is right; `slug`/`title` consistent between
  `meta.yaml` and `models_entry.yaml`; `check_metadata.py` clean.
- **No figure digitised, no review gate outstanding.** Correct.
- **Table 2 transcription** verified against the render: water 998/1/72,
  tetradecane 763/2.2/27, paraffin A 795/2.3/28, paraffin B 790/2.9/28,
  polyacrylamide zero-shear 50/100/190/350. The "(Est.)" heading on the surface
  tension column is carried into the sidecar. Newtonian exclusion is justified.
- **Figure 11 caption** verified: "In the model calculations take Udf = Utrans;
  the experimental data points are plotted with measured Udf values." The page's
  description of the two conventions is exactly right.

---

## VERDICT

**Safe to publish after fixes 1 and 2.** Both are claims the page advertises as
its contribution (they appear in `meta.yaml adds:`/`validation:`, the
`models_entry.yaml` description, the README headline and the notebook's "What
pymrm adds"), and both are wrong as stated: the "complementary closures / sensible
above 14 kg/m3" result is an artefact of a fixed band, and the C1 exponent check
does not check what it says it checks. Findings 3–6 are wording changes in the
same files. 7–10 are optional improvements, though 7 also touches the published
F1.4 page.

The page's headline finding — Eq. 2 collapsing at atmospheric gas density — is
correct, correctly attributed to Krishna and Ellenberger, and independently
confirmed here against three of their own figures and their own notation list.
