# B1.7: The Mears criteria

**Catalog ID:** `B1.7` · **Section:** B · **Tier:** T0 · **Structures:** `S3` + `S4`
· **Provenance tier:** 6 (no measurement, figure or table anywhere in the source)

The experimenter's side of the pellet problem, completing the arc of
`B1.1`/`B1.4`/`B1.6`: the one-line inequalities in **observables** that decide
whether a measured rate is chemistry or transport (intraparticle, interphase,
interparticle), almost all claiming the same promise, *satisfy this and the
rate deviates by less than 5 %*. This page solves the underlying pellet, film
and dispersion problems exactly and measures what each printed boundary
actually admits.

## Provenance

**Mears, D. E. (1971), "Tests for Transport Limitations in Experimental
Catalytic Reactors", *Ind. Eng. Chem. Process Des. Develop.* 10(4), 541-547,
`doi:10.1021/i260040a020`**, on disk at
`~/papers/pymrm-gallery/Mears1971-tests-for-transport-limitations-IECPDD10-541.pdf`,
is the only document consulted.

Two traps the next reader should know:

- **The PDF opens with the tail of the *preceding* article** (a
  catalytic-cracking study discussing "the Voorhies (1945) relationship"; it
  looks like `B2.1`'s paper). The Mears article starts partway down page 1 of
  the file; identity was established from the title block on that page, the
  running head, and the ACS download stamp `article-pdf/10/4/541/`.
- **The text layer destroys exactly the load-bearing symbols**: the
  effectiveness factor extracts as a bare comma, eq. 23's `0.4` as a degree
  sign and `'4`, eq. 25's `5.3` as a degree sign and a dash. Every threshold
  was read off 300 dpi renders, the scan's *native* resolution (CCITT-G4
  bilevel at 300 ppi), so higher renders only interpolate.

**This paper is a review**, and the criteria "derived by the author" are
derived in two companion papers: Mears 1971a (*Chem. Eng. Sci.*, "in press" in
the reference list; the axial criterion) and Mears 1971b (*J. Catal.* **20**,
127; the interphase and interparticle heat criteria). Neither companion was
consulted; the page audits the criteria exactly as the review states them,
which is also how the literature cites them.

## What the page does

1. **Checks the transcription against the paper's own identities**: 0.22 =
   1/sqrt(20), "four and a half times" = sqrt(20), Carberry's 0.1 = 2 x 0.15,
   eq. 15's Biot threshold 10 = 2 x 0.75/0.15, and the perturbation slope n/15
   that makes 0.75 the constant the paper says it rounded to 1.
2. **Audits each criterion at its own printed boundary against exact solves.**
   Eq. 17 is exact at first order (5 % to all digits). The intraparticle family
   is liberal by its own rounding: eq. 1 admits 6.47 % <!-- claim: eq1_dev_at_bound_n1 -->,
   eq. 3 from 6.29 % <!-- claim: eq3_dev_at_bound_n3 --> (n = 3) to
   6.79 % <!-- claim: eq3_dev_at_bound_n05 --> (n = 1/2). Eq. 9 over-certifies
   even with the unrounded constant for gamma*beta in
   (0.2222 <!-- claim: eq9_nopower_window_lo -->,
   1.1586 <!-- claim: eq9_nopower_window_hi -->) at gamma = 30, where
   Petersen's 13 admits 37.8 % <!-- claim: petersen13_dev_gamma30 -->. Eq. 21
   fails structurally at low Biot number, and Mears' axial criterion holds
   everywhere tested: 3.91 % <!-- claim: eq30_dev_X02 --> to
   4.95 % <!-- claim: eq30_dev_X099 --> at first order.
3. **Measures what stacking tolerances costs**: passing eqs. 1, 17 and 21
   simultaneously still admits 11.15 % <!-- claim: joint_pass_worst_dev -->, at
   an interior optimum Bi_m = 6.33 <!-- claim: joint_pass_worst_bim --> where
   two criteria bind at once.
4. **Places `B1.4`'s fold finding**: on the folded branch of gamma*beta = 24,
   eq. 1's band certifies a deviation of
   3135 <!-- claim: eq1_fold_worst_certified_dev -->, while eq. 9's band edge
   admits 7.20 % <!-- claim: eq9_fold_worst_certified_dev --> (a root-found
   supremum), because its bound has shrunk below the ignited branch's reach.

## Validation

No experiment: tier 6. The evidence: three routes to the headline number
(closed form, shooting, pymrm), whose largest gap is the pymrm grid error on
800 cells, 1.04e-7 <!-- claim: headline_three_route_spread -->; pymrm
against the shooting reference at eight threshold states (Dirichlet and film
surfaces, isothermal and Weisz-Hicks kinetics) at second order, worst
8.2e-8 <!-- claim: pymrm_vs_shooting_iso_n800 --> on 800 cells, with the
finest difference at least 1e3 times the round-off of eta; the pymrm
dispersion vessel (van Leer TVD by `ConvFluxTVD` and `newton`) against the
classical closed form, 7.11e-8 <!-- claim: pymrm_axial_vs_closed_form --> on
1600 cells; the two solvers of each discrete system (exact Jacobian and
`NumJac`; Newton and deferred correction) against each other; the paper's own
printed numbers reproduced from other printed numbers; and a break table in
which every claim but three structural identities has its own injected-defect
row that fails it.

## Corrections

- 2026-11-01: said pymrm agrees with the shooting reference at the isothermal
  5 % states to 8.5e-08 on 800 cells at order 1.95. The shooting reference was
  integrated at rtol 1e-10, and the 5 % state root-found on it carried an error
  of 3.6e-9 in eta (reaction order 2), which both inflated the error and bent the
  order. With the reference at rtol 1e-12 the worst error is
  8.2e-8 <!-- claim: pymrm_vs_shooting_iso_n800 --> and the observed orders
  are 2 <!-- claim: pymrm_vs_shooting_iso_order -->, rounded, over all eight states.
- 2026-11-01: said B1.4's eta = 1.0483 at Phi = 1 (inhibited
  Langmuir-Hinshelwood rate) was reproduced by an independently written branch
  tracer. Both pages use the same shooting construction (series start, DOP853,
  stop at y = 1), so the agreement checks the two pages' bookkeeping, not the
  method.
- 2026-11-01: said a final cell maps every one of the 53 agreement metrics to
  the break row that moves it. The cell assigned labels; rows were run on 8 of
  the 53 metrics. Every claim except the three structural identities now has its own row, printed in the break table.

## Files

```
index.py        the page source (jupytext percent format), paired with index.ipynb
index.ipynb     the executed notebook
model.py        rate laws, the shooting reference and threshold finder, the pymrm
                pellet, the film closed forms, the dispersion vessel
test_model.py   fast checks of model.py and two stored claims
meta.yaml       page metadata
agreement.json  the claims, format 2, written by the last cell
```

No `data/`: nothing is digitised, no CSV is shipped, and no other page's
dataset is loaded; three committed `B1.4` claims are read with `load_claims`
and printed beside this page's values.
