# J2.6: Smoluchowski (1917), coagulation kinetics and Zsigmondy's counts

M. v. Smoluchowski, "Versuch einer mathematischen Theorie der Koagulationskinetik
kolloider Lösungen", *Z. physik. Chem.* **92**, 129-168 (1917); the article has no
DOI. The page rebuilds the closed forms of rapid coagulation (eqs. 23, 24, 29),
recomputes Section IV's arithmetic from the printed counts of Zsigmondy's Versuche D,
E and F, and tests the paper's reading that the particles must nearly touch
(a capture radius of the order of the particle diameter; contact is R = 2r) against how precisely the counts fix R/r. A closed-form page;
it calls no pymrm operator.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/smoluchowski-1917-printed-numbers.csv` | the three count tables, the NaCl series and the stated numbers, as printed |
| `data/check_values.csv` | the page's functions at one test condition, computed by a second route |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it.

## Things to know before reusing this

- On the scatter of three rows each, the counts give
  R/r = 1.405 <!-- claim: r2_D_Rr --> ± 0.034 <!-- claim: r2_D_se --> (D),
  3.10 <!-- claim: r2_E_Rr --> ± 0.18 <!-- claim: r2_E_se --> (E) and
  2.617 <!-- claim: r2_F_Rr --> ± 0.053 <!-- claim: r2_F_se --> (F). The choices of the
  conversion (the law inverted, law (16) for (24)_1, the initial count, the rows left out,
  the delay, the kernel, the constant C) move R/r far more: D 1.30 <!-- claim: r6_D_lo --> to
  2.03 <!-- claim: r6_D_hi -->, E 2.27 <!-- claim: r6_E_lo --> to
  4.61 <!-- claim: r6_E_hi -->, F 2.35 <!-- claim: r6_F_lo --> to
  3.58 <!-- claim: r6_F_hi -->. Quote the two together. Laws (11) and (14) are left out
  of these ranges: they hold the nucleus fixed (D for D1 + D2, p. 143).
- Point IV of the paper is a statement of order, and the counts support it: across the
  choices R/r lies between 1.30 <!-- claim: r6_all_lo --> and
  4.61 <!-- claim: r6_all_hi -->. Contact as a sharp value is a separate question: on each choice's
  own standard error, contact is excluded under the paper's conversion in
  3 <!-- claim: r2_n_excluding_2 --> of the three Versuche and in
  19 <!-- claim: r6_n_other_excluding_2 --> of the
  28 <!-- claim: r6_n_other_choices --> other choices. D lies below contact in
  9 <!-- claim: r6_D_n_same_side --> of its 10 <!-- claim: r6_D_n_choices --> choices,
  which the theory does not allow (R/r = 2 is contact, p. 151).
- The paper prints neither the viscosity, the temperature nor Avogadro's number. The
  constant C of eq. 26 implied by its tables is used; water at 20 C gives
  0.986 <!-- claim: r1_C_A1_over_C_P --> times it.
- The printed sum of p. 156 has kt/(t+T) where the sum of eq. 24 gives kT/(t+T), and
  eq. 29 holds with x = t/(kT), not the printed t/T.
- Kernel (18) with compact aggregates moves R/r by -0.43 % <!-- claim: r5_D_shift --> (D)
  to -0.71 % <!-- claim: r5_E_shift --> (E); the transient term
  the paper drops is at most 3.5 % <!-- claim: r4_transient_ratio_max --> of the shift
  the early rows would need.
