# J1.4: IAST 1965, mixture adsorption from pure-component isotherms alone, and a paper that prints no numbers

The ideal adsorbed solution theory is one sentence, the partial pressure of
an adsorbed component is its adsorbed-phase mole fraction times the pressure
it would exert as a pure adsorbate at the same temperature **and spreading
pressure**, $P y_i = P_i^{\circ}(\pi)\,x_i$, and from it, mixture adsorption
follows from pure-component isotherms with no mixture data and no mixture
parameters.

**Source.** Myers, A. L. and Prausnitz, J. M. (1965). *Thermodynamics of
Mixed-Gas Adsorption*. AIChE Journal **11**(1), 121–127,
[doi:10.1002/aic.690110125](https://doi.org/10.1002/aic.690110125). **It is
the only document read for content.** Identity confirmed from its own first
page on a native-resolution render: title, by-line "A. L. MYERS and J. M.
PRAUSNITZ / University of California, Berkeley, California", the
Myers-at-Penn footnote, the running feet, and the Wiley margin naming the
same DOI. `pdfimages -list` reports all seven pages as CCITT-G4 bilevel at
**300 ppi native**; every numeric was read off a digit-scale crop at that
resolution and the text layer was used only as a search index. (PDF page 7
carries the start of the *next* article below the reference list; nothing was
read from it.)

- `index.py`, the page, paired with `index.ipynb` (`jupytext --sync`).
- `model.py`, the isotherms, the IAST solver and the helpers the page calls;
  `test_model.py` checks it offline in seconds.
- `data/myers-1965-printed-claims.csv`, the paper's printed scalars:
  conditions, the text's own counting, reference-list numbers. **The paper
  prints nothing else numeric.**
- `data/iast-illustrative-reference.csv`, **computed by the page** (sidecar
  says so): IAST solutions on two illustrative Langmuir pairs, byte-identical
  across runs, for regression-testing independent IAST implementations.

**No other page's dataset is loaded**, so none of the cross-page
reconciliation obligations apply. **No curve is digitised, no figure is
reproduced, no page image is committed anywhere.**

## The scope decision, made before any code

The paper's four experimental validations, CH₄–C₂H₆ and C₂H₄–CO₂ on
activated carbon, CO–O₂ and C₃H₆–C₃H₈ on silica gel, are **figures**
(Figs. 3–10) whose points belong to four external papers, none on disk. There
is **not one table in the paper** (the only tabular matter is NOTATION and
LITERATURE CITED), no numerical result of the theory is printed anywhere, and
even the pure-component inputs exist only as Fig. 2's *already-integrated*
$\pi A/RT$ curves ("CALCULATED FROM EXPERIMENTAL ADSORPTION ISOTHERMS OF
SZEPESY & ILLES (1963)"). Reproducing Fig. 3 would mean differentiating a
digitised integral of absent data, figure-derived inputs feeding
figure-derived targets, so the experimental case is **scoped out**, and the
page says plainly that it establishes what the theory *is*, not that it is
empirically right. If the Szepesy & Illés tables (Acta Chim. Hung. **35**,
37/53/245) are acquired, the four comparisons become a natural companion
page.

## What the page establishes

1. **The derivation chain closes, with two printed defects in it.** Eleven
   sympy identities, all exactly zero: (11)+(12) → (13); (22)+(41)+(42) →
   (29); (20) → (21) at constant $P$; eq. (19) on Langmuir; the Henry chain
   (30)→(32)→(35) recovering $n_i = K_iPy_i$ (31) exactly; both lever rules;
   the equal-capacity closed form; the extended-Langmuir curl and its
   vanishing at $m_1 = m_2$. Found on the way: **eq. (11) is printed without
   the `+` between its two $RT\ln$ terms** (the product form is dimensionally
   incoherent and provably does not yield eq. (13); the lost `+` is an
   inference, labelled), and p. 123's "With the pressure $P$ of the mixture
   held constant, Equation (21) becomes" **introduces eq. (21) by its own
   number** where the derivation needs (20). Reported, not repaired.
2. **The solver is proved against closed forms it never touches.** The page
   derives symbolically that equal-capacity Langmuir IAST *is* extended
   Langmuir; the Newton solver reproduces it to its own tolerance, below
   1e-12 <!-- claim: equal_m_closed_form_max_rel -->, over a grid of states
   spanning three decades of pressure. In the Henry limit, the case the
   paper itself proves rigorous, the solver is exact on pure Henry isotherms
   (to below 1e-12 <!-- claim: henry_model_solver_max_rel -->) and converges
   to $K_1/K_2$ at first order in $P$: observed orders from successive
   differences at least 0.95 <!-- claim: henry_convergence_order -->
   and at most 1.05 <!-- claim: henry_convergence_order --> on a ladder of pressures a
   decade apart, so the deviation of 2.75e-8 <!-- claim: henry_selectivity_dev_P1em8 -->
   at $P = 10^{-8}$ is a first-order approach to the limit, not a round-off
   figure.
3. **The paper's consistency charge, quantified two ways.** P. 125 reports
   the Langmuir mixture model "was not thermodynamically consistent". Around
   a closed loop in $(P, y_1)$, IAST closes the Gibbs one-form to round-off,
   below 1e-12 <!-- claim: loop_circulation_iast_abs -->;
   unequal-capacity extended Langmuir leaves **11.9 % <!-- claim: loop_deficit_share_of_scale -->
   of the loop's own scale**, by boundary line integral and by the
   sympy-derived curl $(m_2{-}m_1)b_1b_2P^2/D^2$, routes sharing no assembly,
   agreeing to below 1e-12 <!-- claim: circulation_route_rel_agreement -->, and
   exactly linear in $m_2-m_1$ (doubling $m_1 - m_2$ doubles the deficit to
   below 1e-12 <!-- claim: circulation_linearity_dev -->). The shortcut model's
   cost on the illustrative pair: 21.5 % <!-- claim: unequal_m_nt_dev_windowsup -->
   in total amount (window supremum, a bounded maximisation refined from a
   scan; the grid read printed beside it is
   0.67 % <!-- claim: unequal_m_nt_dev_grid_low_frac --> low) and
   64.2 % <!-- claim: unequal_m_n1_dilute_limit_dev --> on the dilute component
   (closed-form limit).
4. **A structural result the paper stops short of** (searched: the Fig. 10
   discussion and the conclusions do not state it): IAST selectivity is
   **exactly constant** whenever the two pure isotherms are the same curve
   shifted in $\log P$, verified on a dual-site pair to the nested solver's
   tolerance (relative spread below 1e-10 <!-- claim: shape_translation_selectivity_rel_spread -->)
   and moved to 1.1e-2 <!-- claim: shape_perturbed_selectivity_rel_spread --> by
   perturbing one site's affinity by five percent, so composition-dependent selectivity
   under IAST is a *shape* signal. The paper's Henry-limit constancy is the
   linear special case.

## Printed defects, all reported, none repaired

- **Eq. (11)**: missing `+` (above), proved from eqs. (9), (10), (12), (13).
- **Eq. (21) cross-reference**: cites itself where (20) is needed.
- **Fig. 6 legend**: "EXPERIMENTAL POINTS OF MARKHAM & BENTON (1950)" [sic]:
  the paper's only Markham & Benton reference (ref. 8) is *J. Am. Chem. Soc.*
  **53**, 497 (**1931**), and the p. 125 text credits the 0 °C data to "(8)";
  1950 is ref. 7 (Lewis et al.).
- **Ref. 5**: "Ibid, 63, 456 (1959)", the Ibid (= *J. Chem. Phys.*) cannot
  carry vol. 63 in 1959 against ref. 4's own (17, 1949): 4.6 <!-- claim: ref5_volumes_per_year --> volumes/year.
  *J. Phys. Chem.* **63** is 1959; that it was intended is an inference,
  labelled. Ref. 5 is cited nowhere in the running text (every "(5)" is
  equation (5); text-layer search, hits verified on crops).
- **"The absorbent is not counted as a component"** [sic, p. 122]:
  "adsorbent" is spelled correctly everywhere else.
- **Serpenskii / SERPINSKII**, text and reference list vs the three figure
  legends: one surname, two transliterations (an observation, not an error
  claim).

## What pymrm adds

Nothing to the theory. The whole $(P, y_1)$ grid is solved as one Newton
problem (`newton` on `NumJac((K, 1))`, iterating in $\ln z$ to stay off the
$P^{\circ}(0) = 0$ singularity), and every solve raises unless `newton`
reports success; the dual-site pair's isotherm inverse is a nested vectorised
`newton`, the machinery tabulated-isotherm users actually need;
`pymrm.checks.observed_orders` gives the orders of the two quadratures and of
the Henry-limit approach from successive differences, each with its round-off
margin. Extrema are closed forms or bounded maximisations refined from a scan,
never sampled. The most reusable output is not pymrm at all: the
**closed-loop Gibbs consistency diagnostic**, applicable to any proposed
mixture model $n_i(P, y)$ with no IAST solve required.

## Corrections

- 2026-11-17. The page said the shortcut model's total-amount error was a
  "21.5 % in total amount (window supremum, root-found; the grid read printed
  beside it is 0.67 % low)". The value stands
  (21.5 % <!-- claim: unequal_m_nt_dev_windowsup -->), but it was never
  root-found: the maximum over $y_1$ came from a bounded Brent maximisation.
  It is now a bounded maximisation refined from a scan of the whole $y_1$
  range, raising when the scan maximum sits on a scan end.

## Caveats worth reading before reusing anything

- **Nothing here is experimental and nothing is a fit.** Tier 6 by
  construction: pairs L, U, S are illustrative parameter sets in arbitrary
  units. The *magnitudes* (11.9 % <!-- claim: loop_deficit_share_of_scale -->,
  21.5 % <!-- claim: unequal_m_nt_dev_windowsup -->,
  64.2 % <!-- claim: unequal_m_n1_dilute_limit_dev -->, factor
  4.39 <!-- claim: selectivity_dilute_ratio_P1 -->) belong to
  those parameters; the *structural* statements (deficit nonzero iff
  capacities differ, exactly linear in $m_2-m_1$; selectivity constant iff
  shapes match) are the claims.
- **Eq. (14)** (fugacity form) is transcribed and deliberately not exercised:
  the paper prints no real-gas model to put in it.
- **The four data sources are not consulted**; the claim that Arnold's model
  is inconsistent is the paper's, quoted as such. The inconsistency this page
  *proves* is extended Langmuir's, from its own closed form.
- **The equal-capacity closed form and the shape-translation result are this
  page's derivations** (standard results; no source was consulted for them),
  proved symbolically on the page.
- **This page is not a rate.** Breakthrough is `J1.5`; this page's `iast`
  drops into that page's mixture-equilibrium slot. `J1.1` (Langmuir) and
  `J1.3` (BET) supply pure isotherms for the $F_i$ slot; both were read,
  neither is loaded, and no number of theirs is retyped here.
