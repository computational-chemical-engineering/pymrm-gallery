# A1.8 — Two ways to close gas–solid drag

One scalar decides almost everything a two-fluid fluidised-bed simulation does,
and the MFIX Theory Guide says in one paragraph why there is no single answer to
it. This page builds the two families that paragraph names, divides them, and
puts each against a measurement it never saw.

- **Structures:** `S3` — algebraic closures. **No pymrm operator appears**, and
  the page says so where a reader would look for a solver, as `A1.1` and `A1.6`
  do. (The catalogue says `S1`; changed deliberately to match the section-A
  siblings. `check_metadata.py` will not catch this — it compares `meta.yaml`
  against `models.yaml` only.)
- **Runtime:** ~4 s
- **Data tier:** 2 — two tables printed in papers, neither closure fitted to
  either.

## Source

**Syamlal, M., Rogers, W. and O'Brien, T. J.**, *MFIX Documentation Theory
Guide*, Technical Note DOE/METC-94/1004 (DE94000087), Morgantown Energy
Technology Center, U.S. Department of Energy (December 1993),
[doi:10.2172/10145548](https://doi.org/10.2172/10145548) — 54 pp, freely
downloadable from OSTI. Section 2.2.1, journal pages 10–11, read on **400 dpi**
renders, which is the embedded bilevel page images' own native resolution.

**Do not use this file's text layer for any digit.** It renders `0.06 Re` as
`0.O6Re` with a capital letter O, ε_g^2.65 as `_g-2'6s`, "Ergun (1952)" as
`Ergun (f952)`, and the equation numbers (12) and (16) as `(121` and `(1''`.

**The origin of eq. (11) is an unpublished report.** The report's own reference
list gives Syamlal and O'Brien (1987), *"A Generalized Drag Correlation for
Multiparticle Systems," Unpublished report* — so this DOE report is not a
convenient source for the terminal-velocity-to-drag conversion, it is the citable
**published** one. The catalogue's citation for the case, "Syamlal & O'Brien
(1989)", is a *different* item in the same list: "Computer Simulation of Bubbles
in a Fluidized Bed", AIChE Symposium Series No. 270, **85**, 22–31, a
bubble-simulation paper. Both were read off a 400 dpi render of report page 44.

**Origins cited but not consulted.** Garside & Al-Dibouni (1977) for eqs.
(12)–(14) and Dalla Valle (1948) for eq. (16); both are printed in full with
attribution in the 1993 report, which is the reprint route of `AGENTS.md`.

**Also read directly:** Richardson, J. F. and Zaki, W. N., *Sedimentation and
fluidisation: Part I*, Trans. Instn Chem. Engrs **32** (1954) 35–53, from the
verbatim Golden Jubilee reprint, Trans IChemE **75** (Dec 1997) S82–S100,
[doi:10.1016/S0263-8762(97)80006-8](https://doi.org/10.1016/S0263-8762(97)80006-8),
on 300 dpi renders. Its Table VI is shipped here. That paper is `A1.5`'s source;
`A1.5` is untouched and only one of its tables is used.

## Scope — read this before anything else

The catalogue asks for **Gidaspow / Syamlal–O'Brien / Wen–Yu**. The document on
disk carries **Syamlal–O'Brien complete and neither of the other two**:

- the strings "Wen and Yu" and "Wen & Yu" do **not occur anywhere** in the
  report;
- the report prints **no Gidaspow drag closure and no blend rule**. Gidaspow
  (1986) — "Hydrodynamics of Fluidization and Heat Transfer: Supercomputer
  Modeling", Appl. Mech. Rev. **39**, 1–23 — is cited exactly once in §2.2.1, for
  the observation that a packed-bed correlation must be supplemented at low
  solids fraction.

  **Do not restate that as "Gidaspow is cited only in passing" — an earlier
  version of this page did, in four files, and it is false.** Measured on the
  extracted text layer, the name occurs **28 times** (18 in the body, 10 in the
  reference list), across four Gidaspow-first-author entries — 1974 Round Table,
  Gidaspow & Ettehadieh 1983, Gidaspow et al. 1984, Gidaspow 1986 — plus Ding &
  Gidaspow 1990, Syamlal & Gidaspow 1985, Shi/Gidaspow/Wasan 1987,
  Arastoopour/Lin/Gidaspow 1980, two Lyczkowski/Gidaspow papers and Bouillard et
  al. 1989. More than citation: the report **adopts a Ding & Gidaspow (1990)
  expression as its own eq. (88)** (granular-energy transfer, §2.5.3) and a
  Syamlal & Gidaspow (1985) conductivity model in the heat-transfer section. The
  text layer mangles names, so that count is a floor, not a ceiling. What is
  load-bearing, and true, is only that no Gidaspow *drag* law and no blend rule
  appear anywhere.

Nothing about either is written here from memory. What replaces them is the pair
the report's own paragraph on journal page 10 points at, both on disk: **Ergun
(1952)**, via `A1.1`, and **Richardson & Zaki (1954)**. Two documents would
complete the case, and they are named on the page and in `models_entry.yaml`:
Gidaspow's 1994 monograph, and the Wen & Yu Chem. Eng. Progr. Symposium Series
paper (**not** the A.I.Ch.E. Journal communication on disk, which is `A1.6`'s
source and contains no drag law).

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 4 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, never the notebook |
| `data/mfix-1993-syamlal-obrien-constants.csv` | the **eleven** printed constants of eqs. (11)–(16), each tagged with its equation number |
| `data/richardson-zaki-1954-table6.csv` | Table VI — **18** tabulated velocity–voidage indices over Re 0.39–489, plus the paper's own two log₁₀ columns. Each n₀ is a fitted slope, then extrapolated along a second fitted line, so "tabulated", not "measured" |
| `meta.yaml` | page metadata |

Both CSVs have a `.meta.yaml` provenance sidecar. Three further datasets are
**borrowed**: `A1.1`'s Ergun constants and 244 markers, and `A1.7`'s Geldart
table. The page's *The data* section lists what those pages have already
established about the rows it uses.

## What the page checks, in the brief's order

| route | check | result |
|---|---|---|
| 2, identity | `V_rm(eps_g = 1) = 1` at every Re | **2.2e-16** — the radicand of eq. (12) is a perfect square. **Blind to all three voidage exponents, to the 0.8 and to the switch voidage** (1^x = 1); 5 of 10 injected mis-readings leave it at the threshold |
| 2, identity | the force balance returns `v_slip = V_rm v_t` | **9.3e-15** over four decades of size. Catches the `V_rm^2` (23×), the `C_Ds` argument (98 %), the 3/4 (28 %) and the `eps_g` (63 %); **blind to every constant of eqs. (12) and (16)** |
| 3, measurement | vs Richardson & Zaki's Table VI, 18 indices, nothing fitted | implied index biased **+3.2 %** at eps_g = 0.45, **+5.4 %** at 0.50, **+29.2 %** at 0.95 — grows with voidage because the closed form is not a power law. Viscous limit is a one-liner: 4.14 + 1 = **5.14 vs their tabulated 4.65, +10.5 %**. **Conditional on a reading — see below** |
| 3, measurement | vs 8 minimum-fluidisation velocities Geldart measured in 1973 | **all three overpredict at the derived voidage** — Ergun **+42.1 %**, Syamlal–O'Brien **+101.2 %**, Richardson–Zaki **+111.1 %**. No closure is validated; the *ordering* is stable and Ergun is least biased at every voidage tried |
| 2, identity | Table VI against its own two log₁₀ columns | **8.0e-04** and **5.4e-04**; a single-digit slip in either moves it to 1.00 |
| 2, cross-page | `A1.6`'s and `A1.1`'s published numbers, recomputed by independent code | Ergun at eps_MB **+58.43 %** vs its 58.435; at the settled eps₀ **+42.06 %** vs its 42.060; at a textbook 0.42 **+4.82 %** vs its 4.8185; Wen & Yu eq. (1) **−25.21 %** vs its −25.214; `A1.1`'s refit reproduced to **2.6e-06 %** (the worse of the two constants; k₁ alone is 4.6e-07 %) |
| 4, digitised | — | **nothing on this page is digitised.** The one digitised dataset it loads is `A1.1`'s, used only to bound Ergun's own constants |

**The Check-3 break table is swept over voidage, and that is a correction.** It
used to be run at eps_g = 0.50 only, and the page claimed on that basis that
"between the two checks every printed constant of the closure is exercised by
something". **That claim was false.** Eq. (14) is branched, so the evaluation
voidage decides which of its constants are live: at 0.50 the dilute branch is
never reached, and `B_exponent_dilute` (2.65) and `B_switch_voidage` (0.85) moved
*nothing* — in any of the four checks. Check 1 sits at eps_g = 1 where 1^x = 1,
Check 2 never reads eq. (14), and Check 4's voidages are 0.441–0.498, all dense.
Meanwhile the page's own headline, +29.2 % at eps_g = 0.95, is computed **on** the
untested branch: mis-reading 2.65 as 1.65 moves it −15.3 points. The check that
would have caught the error was the one the page published as its result. Swept
over six voidages, 7 of 15 injected mis-readings that the 0.50 column cannot see
become visible. This is the third time in this repository that a page's guard
structure has contained the defect it exists to catch, and the page says so.

## Headline results

| | |
| --- | --- |
| the two closures at eps_g = 0.44, viscous limit | **1.192** — they agree to 19 % |
| the same, inertial limit | **0.957** — 4 % |
| the same, worst point in between | **0.579 at Re_m = 4.09**, i.e. a factor **1.73** — but on Ergun's own abscissa that is Re/(1−ε) = **7.30**, and his 244 recovered markers start at **7.45**, so the worst point sits 2.0 % *below* his data. At eps_g = 0.50 the dip is at 9.72, comfortably inside |
| where the dip stops being an interior minimum | above eps_g ≈ **0.50**; above that they simply diverge |
| Ergun drag / isolated-sphere drag as eps_s → 0 | **7.9 %** at Re = 1, **0.10 %** at Re = 0.01, **5.79×** at high Re |
| the eps_g = 0.85 branch switch in eq. (14) | **+0.049 %** in B; at most **0.049 %** in V_rm. **It is not a discontinuity** |
| Dalla Valle's viscous limit | C_D·Re → **23.04** against Stokes' 24, **−4.00 %** |
| implied index vs Richardson & Zaki, eps_g = 0.45 → 0.95 | **+3.2 % → +29.2 %**, on a 15-point row-to-row spread — **under the slip reading of V_rm**. Under the superficial reading it is **−26.3 % → −0.3 %** |
| Geldart's 8 spherical cuts, at the derived eps₀ | Ergun **+42.1 %**, Syamlal–O'Brien **+101.2 %**, Richardson–Zaki **+111.1 %** |
| the same, at a textbook eps_mf = 0.42 (`A1.6`'s choice) | Ergun **+4.8 %**, Syamlal–O'Brien **+59.6 %**, Richardson–Zaki **+56.8 %** |
| voidage each closure would need to be unbiased on those rows | Ergun **0.4147**, Richardson–Zaki 0.3708, Syamlal–O'Brien 0.3538 — the last two below random loose packing |
| Syamlal–O'Brien / Ergun on those rows | **0.97 → 1.68**, rising over the first seven cuts and **turning over on the eighth** (1.6812 at 263 µm, 1.6717 at 318 µm). **Not** monotone |
| Syamlal–O'Brien / Richardson–Zaki on those rows | **0.74 – 1.12** — an *algebraic identity*, not a measured result: v_t cancels out of it to 2.2e-16 |
| Wen & Yu's empirical u_mf correlation, same rows | **−25.2 %** — and it is the one comparison with no voidage in it at all |

## What is not here

- **Nothing about Gidaspow's closure or Wen–Yu drag.** The eps_g = 0.8
  discontinuity a comparison of blended closures exists for is neither computed
  nor asserted. The only switch measured here is Garside & Al-Dibouni's at 0.85,
  and it turns out to be continuous.
- **No pymrm operator.** Both closures are algebraic and the two measured
  comparisons are roots of scalar equations. Inventing a PDE would obscure the
  comparison.
- **The inertial half of every closure is essentially untested against
  measurement.** Every Geldart row is a group A powder at Re_t < 28 and a bed
  Reynolds number below 2. Table VI reaches Re = 489 but constrains eqs.
  (12)–(14) only — Check 3 is exactly blind to the single-sphere drag curve.
- **No claim that Ergun is a reference standard, and no claim that any closure is
  validated.** `A1.6` shows Ergun is +58 % biased at Geldart's reported ε_MB and
  +4.8 % at a textbook 0.42; both are recomputed here rather than quoted. What
  the page *does* now support, after an earlier version refused to rank at all,
  is the narrow statement that **Ergun is the least biased of the three at every
  voidage tried** — with the absolute bias attributed to the voidage, which
  nobody measured, rather than to the closures.
- **Neither Check-4 ratio is measurement-facing.** v_t cancels identically out of
  Syamlal–O'Brien / Richardson–Zaki (verified to 2.2e-16), so 0.74–1.12 is
  Check 3 restated at Geldart's eight (ε₀, Re_t) pairs and would read the same if
  the experiment had never been done. Geldart's U₀ enters only the per-cent
  columns.
- **Which of two readings of V_rm the report meant.** Eq. (11) forces the slip
  ratio; the report's own prose introducing V_rm implies the superficial one.
  They are one unit of index apart and the alternative **reverses the sign** of
  Check 3 (−26.3 / −24.2 / −0.3 %, viscous limit −11.0 %). Both are printed and
  CI-tracked; the page adopts the one eq. (11) implements and says so. Garside &
  Al-Dibouni (1977), which would settle it, is not on disk.
- **n_table6 holds Table VI's index constant outside Re 0.39–489.** Four of the
  21 context rows fall below the bottom; none of the eight Diakon rows does.
- **Table VI is not a raw measurement.** Each n₀ is a fitted slope, extrapolated
  along a second fitted line to d/D = 0, from liquid–solid runs. And Garside &
  Al-Dibouni fitted their relation to overlapping literature, so the comparison
  measures the cost of one closed form standing in for another, not two
  independent witnesses.
- **The switch voidage 0.85 is the one printed constant no data comparison
  settles** when it is displaced by a little. The Check-3 break table is swept
  over voidage and catches it when the mis-read switch crosses an evaluation
  voidage (+17.5 points at ε = 0.60 for 0.85 → 0.55, +44.5 at 0.95 for
  0.85 → 0.95), but a small displacement is bounded by the +0.049 % continuity of
  eq. (14) and cannot be worth more.
- **The 13 non-Diakon Geldart rows** carry no printed sphericity, and
  Syamlal–O'Brien has no shape factor to give them. Whole-table numbers are
  printed as context, not as a result.

## Regenerating

```bash
python build_page.py                       # rewrite index.ipynb
python - <<'EOF'
import nbformat; from nbclient import NotebookClient
nb=nbformat.read("index.ipynb",as_version=4)
NotebookClient(nb,timeout=1800,kernel_name="python3",resources={"metadata":{"path":"."}}).execute()
nbformat.write(nb,"index.ipynb"); print("OK")
EOF
```

## Related

[`A1.1`](../A1.1-ergun-pressure-drop/) (the packed-bed branch and every Ergun
constant used here), [`A1.6`](../A1.6-wen-yu-minimum-fluidisation/) (the
empirical u_mf correlation, and the caveats that bound Check 4),
[`A1.7`](../A1.7-geldart-classification/) (the measured dataset and the derived
voidage), `A1.5` (Richardson–Zaki in its own right), `E1.2` and `E2.1`
(fluidised-bed models that consume u_mf rather than a drag law).
