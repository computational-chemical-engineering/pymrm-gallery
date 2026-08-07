# J1.1 — Langmuir 1918: the isotherm that could not have proved its own conclusion

Everybody meets this paper as one equation and the straight line that fits it.
Langmuir got the two constants by **drawing a line through the $p/q$-versus-$p$
plot with a straightedge**, so the $q_\mathrm{cal}$ columns of Tables II–XVII are
a two-parameter fit compared against its own training data. That is a goodness of
fit. It is labelled one everywhere on this page, it is never reported as an
agreement metric, and a null baseline is printed beside every agreement quoted
anywhere.

The paper's real evidence is elsewhere, and the page is organised around finding
it.

**Source.** Langmuir, I. (1918). *The Adsorption of Gases on Plane Surfaces of
Glass, Mica and Platinum*. J. Am. Chem. Soc. **40**(9), 1361–1403,
[doi:10.1021/ja02242a004](https://doi.org/10.1021/ja02242a004). Received June 25,
1918. **It is the only document read for content.** Identity confirmed from its
own title page on a native-resolution render: the contribution line
"[Contribution from the Research Laboratory of the General Electric Co.]", the
title, the by-line "By Irving Langmuir", "Received June 25, 1918", and the ACS
download stamp `jacsat/article-pdf/40/9/1361/` printed down the margin of every
page. **PDF page 1 opens with the numbered summary of the preceding article** on
the absorption spectra of metals in liquid ammonia; Langmuir's title sits below
it on the same page.

`pdfimages -list` reports every page as **CCITT-G4 bilevel at 300 ppi native**,
so pages were rendered at 300 ppi and every numeric cropped and re-read at digit
scale. **The text layer was not used for any digit**, and that discipline was
load-bearing: it renders eq. (37)'s coefficient as $25.2\times10^{16}$ where the
page prints $25.2\times10^{\mathbf{15}}$ — a factor ten in every surface coverage
on the page, which would put **17 of 20** entries above the monolayer bound and
reverse the paper's conclusion. (That count is not typed anywhere: the notebook
counts it from the data in §6.5 and the break table prints the same 17.) It also mangled half the $b$ column of Table X.
One ambiguous digit (Table VI, $p = 4.6?$) was settled by arithmetic, not by
pixel shape: only $4.65$ reproduces the printed $q_\mathrm{cal} = 5.45$.

- `build_page.py` — regenerates `index.ipynb` (run from this directory).
- `index.ipynb` — the page. Runtime ~13 s.
- `data/langmuir-1918-isotherms.csv` — every $(p, q_\mathrm{obs}, q_\mathrm{cal})$ of Tables II–XI and XIII–XVII, both bulbs.
- `data/langmuir-1918-constants.csv` — the fitted $a$, $b$ per table joined to the derived columns of Tables XVIII–XX.
- `data/langmuir-1918-lives.csv` — Tables XXI and XXII.
- `data/langmuir-1918-platinum-N0.csv` — Table XXV.
- `data/langmuir-1918-printed-claims.csv` — the 36 scalars living in prose, equations or table headers.

**No other page's dataset is loaded**, so none of the cross-page reconciliation
obligations apply. **No curve is digitised and no figure is used for anything** —
the paper's only figure with content is Fig. 1, a schematic checkerboard of
surface atoms with no numbers on it. **No page image is committed anywhere.**

## Scope, decided before transcription

Forty-three journal pages, twenty-five tables. Most of it is an experimental
campaign, and from journal page 1393 on it is irreversible chemisorption on
platinum with no isotherm in it at all.

This page is about **the isotherm and the constants derived from it**: the theory
of journal pages 1368–1376 (Cases I–VI, eqs. 1–33) and Tables II–XI, XIII–XXII.
Table XXV is carried for one reason — it exercises eq. (37) on a *third*
adsorbent — and Tables I, XII, XXIII and XXIV are **scoped out and not
transcribed**. The decision was taken the way the brief asks: *prefer whatever
lets you test rather than restate*. Restating $q_\mathrm{cal}$ would be restating
a fit.

## What the page establishes

**1. Langmuir stages a fair contest and wins it — 3.82×.** Table VII carries two
extra columns in which he fits Freundlich's $q_F = 8.4\,p^{0.417}$ to the *same*
eleven points with the *same* number of parameters and prints the residuals
beside his own. That is a discrimination, because it can fail. Reproduced with
**both** forms refitted optimally — neither of his two fits is the least-squares
optimum, and a rival should be beaten at its best — the best two-parameter
Freundlich is **3.8235 times worse in RMS** on Table VII and worse on **13 of the
14** tables with four or more observations (median 2.1211). The one exception is
Table XIII, which is exactly the table journal page 1388 says "do not give a
straight line when $p/q_\mathrm{obs}$ is plotted against $p$". Null baselines
beside every ratio: a constant is 19.4× worse than Langmuir on Table VII and a
proportional law 26.8×, and never better than 3.7× and 1.8× anywhere.

**2. The monomolecular conclusion rests on $\beta$, and the margin is 17 %.**
Nothing in the chain $b \to b' \to N_0 \to \beta$ is adjustable: $b$ from the
isotherm, $s$ from a ruler, eq. (37) from Avogadro's number, and the denominator
from **liquid densities**, which contain no adsorption in them. It comes out
below one on all twenty entries, as journal page 1391 claims. The honest margin
is the tightest row: $\beta = 0.857$, i.e. **16.7 % of headroom**. The correction
goes the right way — the paper says the mica figures are 10–30 % too high because
the blank was never subtracted, which puts that worst $\beta$ between 0.600 and
0.771 — and no correction is applied anywhere on the page.

**3. Langmuir's Case VI *is* the BET isotherm, twenty years early.** This is the
page's own result and it is not in the paper. Eqs. (26), (29) and (30), with
Langmuir's own stated assumption $\sigma_3=\sigma_4=\dots=\sigma_2$, become

$$\frac{N\eta}{N_0} \;=\; \frac{c\,x}{(1-x)\bigl(1+(c-1)x\bigr)},
\qquad x = \sigma_2\mu \;(= p/p_0 \text{ by eq. } 39),\qquad c = \sigma_1/\sigma_2,$$

with a `sympy` residual of identically zero. **Five** separate symbolic checks,
each of which can fail and each with a break row that moves it: eq. (26) really
follows from eqs. (22)–(25); eq. (30)'s $a$, $b$ **and** $c$ really are the
coefficients of eq. (29) — three printed expressions the paper never derives;
both of the paper's collapse claims about eq. (30) hold exactly; and the same
occupancy ladder summed **to infinity in closed form**, with no truncation
anywhere, gives $cx/[(1-x)(1+(c-1)x)]$ exactly — so the identity does not depend
on stopping eq. (29) at $b$. A numerical ladder summed term by term, sharing no
algebra with either closed form, returns the same answer to $5.6\times10^{-16}$.

**The collapse is Langmuir's own, and so is its consequence.** Journal page 1375,
printed immediately under eq. (30) and read here on a 300 ppi crop: *"If
$\sigma_1$ and $\sigma_2$ are different, but all subsequent values of $\sigma$
(i. e., $\sigma_3$, $\sigma_4$, etc.) are equal to $\sigma_2$, then all the
coefficients in (29) after b are zero."* He states the assumption **and** that
under it the series terminates — which is what licenses truncating eq. (29) at
$b$. (Journal page 1374 carries the weaker physical motivation for the same
collapse; the page quotes both and rests on the stronger.)

That the two theories are related is **not** this page's discovery and the 1938
paper says so itself: its Section II is headed *"Generalization of Langmuir's
Theory to Multimolecular Adsorption"*, read here on a 300 ppi crop of journal
page 311. What is new is that Langmuir's **own** Case VI closed form already *is*
that generalisation.

**4. And therefore the isotherm could not have decided the question.** Eq. (39)
gives $\sigma\mu = 1$ at saturation, so Table XXI's liquefied-gas lives invert
through eq. (3) into the saturation pressures the paper never prints
($p_0$ for nitrogen at 90 K comes out **3.775 atm**). **Every pressure in the
entire paper lies below $p/p_0 = 1.5\times10^{-3}$**, where Case VI and Case I
differ by at most **0.16 %**. Measured against the observational scatter: the
smallest ratio of scatter to theory separation is **12.2**, and on
Table VII, the best nitrogen set, it is **1104**. The two theories are not
distinguishable by these data by one to three orders of magnitude.

*Two scopes, kept apart on the page.* The separation column needs a fitted
isotherm, so it exists only for the **18** eq. (31) entries of Bulb A — Tables
XIII and XVI use eq. (33) and the blank bulb $A'$ ran higher than $A$ on six
tables. $p/p_0$ needs no fit and no bulb, so it is computed over **every pressure
printed in the paper**: 30 table/gas/bulb entries, 124 pressures, both isotherm
forms. Widening the scope does not move the maximum — the eq. (33) tables reach
$5.1\times10^{-4}$ and the blank bulb the same, against the headline
$1.5\times10^{-3}$ from Table XI.

**And the sharper version of the same test: give Case VI both parameters and let
it try.** Refitting *both* models to Table VII from scratch, two free parameters
each, on a 4000-point log grid with a Brent refinement, the best Case VI (BET)
reaches RMS 0.3501824 against the best Case I's 0.3501868 — the multilayer model,
given complete freedom, improves the fit by **1.26 parts in $10^5$**. And it does
it by *becoming* Case I: the fitted $c = 605823$ against $ap_0 = 605812$, the
fitted $v_m = 38.8932$ against $b = 38.8936$. That confirms the conclusion by a
route that never mentions $p/p_0$ separations, and the break row shows the test
has power — fed a curve Case VI generated at $p/p_0 = 0.3$, the ratio collapses
to $9\times10^{-12}$. The Case I side of it is also a second, independent route
to §6.6's fit (different objective assembly, different search); the two agree to
$2\times10^{-15}$.

**The input a sceptic should attack first is attacked on the page.** The
recovered $p_0$ all come from Table XXI, whose vapour pressures Langmuir does not
print, and the 155 K entries invert to tens of atmospheres for the permanent
gases. The conclusion does not depend on them, and the page measures that rather
than asserting it: the headline maximum comes from Table XI, carbon dioxide at
155 K, whose recovered $p_0$ is the **smallest** in the whole set (0.113 atm) and
whose Table XXI life (16.5 s) is the largest; **discarding every 155 K row
outright** moves the maximum $p/p_0$ *down* to $7.7\times10^{-4}$, the largest
departure down to 0.082 % and the smallest scatter-to-separation ratio *up* to
31.1. And the four 155 K permanent-gas rows — precisely those with the largest
recovered $p_0$ — carry the *smallest* relative pressures in the paper
($1.2\times10^{-6}$ at most), so dropping them can only lower the maximum.

*And the direction that would weaken it, which a sensitivity that only drops rows
never probes.* At 155 K the "liquefied gas" of eq. (39) is a **supercooled**
liquid for CO₂; eq. (39) and BET's $p_0$ both mean the liquid, so the page's
reading is the self-consistent one, but the honest question is how wrong it would
have to be. Root-found on the exact expression: parity (scatter = separation) on
the tightest entry needs the recovered $p_0$ to be a factor **12.1 too high**. At
these $x$ the separation is very nearly proportional to $1/p_0$, so that factor
lands within 0.4 % of the scatter/separation ratio itself — it is **not**
independent information, it is that ratio restated in the units a sceptic of
Table XXI would use, and the page says so.

The ratio of the two isotherms is strictly increasing in $x$ (proved: the
derivative's numerator is $c(c-1)x^2+2(c-1)x+2>0$), so the maximum over a table
is *at* its highest pressure — an argument, not a sampled maximum. One per cent
separation is root-found at $p/p_0 = 0.00989933$ against the closed form
$1/101$ valid as $c\to\infty$; that is **1114 times** the highest pressure
Langmuir reached. Reading it off a 40-point log grid instead is 6.3 % out, and
both numbers are printed in the same sentence.

## Four printed defects, reported and none repaired

**Eq. (15).** Printed as $t = \frac{N_0}{N\nu_1(1+\sigma_1\mu)}
\ln\frac{\theta'}{\theta_1-\theta'}$. Put $\theta'=\tfrac12\theta_1$ into it —
which is exactly what Langmuir says he does two lines below — and the logarithm
is $\ln 1 = 0$, so it gives $t_{1/2}=0$ and *not* eq. (16). It also gives
$t=-\infty$ at $\theta'=0$, where the integration starts. Pinning what is not
free (eq. 14 is an ODE; eqs. 7, 10 and 16 are printed), the only free thing left
is the numerator inside the logarithm, and putting $\theta_1$ there — the loss of
a subscript — makes $t(0)=0$ and returns eq. (16) exactly. The page does not
argue this: it **marches eq. (14)** with a `newton`-solved implicit step and
compares against both branches, then root-finds the half-coverage time **on the
marched trajectory** — 0.110160 against eq. (16)'s 0.109954, a discretisation gap
that halves with $\Delta t$ (1.51e-2, 7.52e-3, 3.75e-3, 1.88e-3 at 200→1600
steps, observed order 1.00). That is the independent confirmation. The cell also
inverts the *closed form* at half coverage and gets $2\times10^{-16}$ — and
labels it an **algebraic identity that could not have failed**, evidence for
nothing, kept only because it closes the algebra of eqs. (7), (10), (14) and
(16). The correction is labelled an inference. Nothing the paper concludes is
affected: it never measures a rate.

**Table IV's $b = 58.3$.** Tables IV, V and XVII have exactly **two**
observations, so the line is not fitted but *determined*. Table IV's two points
solve to $b = 57.4937$ and the paper's own $b' = 100.0$ requires $b = 57.5000$ —
two independent routes agreeing to 0.011 % — while the printed 58.3 is 1.40 %
from both. Tables V and XVII pass the same test. Table IV is also the only mica
table with no $q_\mathrm{cal}$ column printed.

**Table XX's $\beta = 0.36$** for methane on glass, where the paper's own
$N_0 = 0.288\times10^{15}$ and its own methane monolayer count $0.63\times10^{15}$
give 0.4571. Proved twice: the $\beta$ identity is exact on 22 of the 23 rows,
and journal page 1391's claim that the adsorption order "is the same in all three
sets of experiments … as indicated by $b'$, $N_0$ or $\beta$" holds under $b'$
and $N_0$ in all three sets and under the recomputed $\beta$ in all three, but
**breaks in one set under the printed $\beta$**.

**"(100 bars) … 0.20"** on journal page 1384. For eq. (31),
$d\ln q/d\ln p = 1/(1+ap)$ exactly; with Table VII's own $a = 0.156$ the printed
0.684 at 3 bars is reproduced (0.6812), but 0.20 corresponds to **25.64 bars**,
not to the printed 100, where the same curve gives 0.0602 — the printed value is
3.32× larger. Two readings are printed and neither adopted.

## And, kept separate from those, one reporting convention

**The $\sigma$ column** of Tables XVIII–XX runs systematically low: eq. (38)
reproduces all 21 cells to 3.77 %, and **18 of the 21** printed values lie below
the computed one, mean $-0.77$ %. That pattern is real and is not a coin toss —
but it is **not an error against Langmuir**, and the page files it separately for
that reason. The three worst rows (−3.77, −2.40, −1.91 %) are exactly the rows he
prints to *one or two* significant figures (10,000; 10,000; 26,000), and
truncating the computed value at each cell's own printed granularity reproduces
15 of the 21 cells against 11 for rounding. He discarded figures; he did not
miscompute them.

## What pymrm is doing here, honestly

Nothing to the isotherm. Eqs. (9), (26), (29), (31) and (33) are closed forms in
one variable; there is no grid, no time step in space, no boundary condition and
no transport, and most of the notebook would run with pymrm uninstalled —
exactly as on `A1.6`, `A1.1` and `J1.3`. `newton` and `NumJac` earn three narrow
places, each because the cheap alternative is *wrong* rather than slow: a fit
that shares no algebra with Langmuir's linearisation — that independence is in
the *objective* (on $q$, not on $p/q$), and §7.2 carries the break row that
refits his way and moves the discrimination metric. The amplitude is profiled out
and the remaining stationarity condition is root-found twice, with `newton` and
with Brent, agreeing to $3\times10^{-11}$ over 29 fits; **that** number is a
root-finder cross-check and nothing more, because both solvers are handed the
same objective, the same profiled amplitude and the same central-difference
derivative and differ only in how they iterate. The page says so where it prints
it. Then eq. (14) integrated rather than quoted, which is how eq. (15) is
disproved; and root-finds where a sweep would have been wrong. The fourth thing is not pymrm at all and is the most useful
output on the page: the Case VI ≡ BET identity, which needed `sympy` and a change
of variables the paper itself supplies.

The page deliberately does **not** fit eq. (29) to Langmuir's data. At
$p/p_0 < 1.5\times10^{-3}$ the parameter $c$ is not identifiable and any fit
would return the Case I answer with an arbitrary $c$. Saying so is the finding.

## Caveats worth reading before reusing anything

- **Tier 5, and only one column of it is measurement.** $q_\mathrm{obs}$ is
  measured; $a$, $b$, $b'$, $\sigma$, $N_0$, $\beta$ and Tables XXI, XXII, XXV
  are all constants Langmuir *derived* from those measurements. Reproducing one
  is reproduction, never validation.
- **The mica numbers are uncorrected** (10–30 % too high, the paper's own words),
  and one Table IX point is disowned by the paper; it is kept, flagged, and the
  discrimination reported with and without it.
- **One external input, and it is what makes $\beta$ a test.** The monolayer
  counts come from liquid densities, and only nitrogen's molecular volume (35.5
  cu. cm.) is printed. The other four are **recovered** by inverting the printed
  counts and labelled recoveries; nothing depends on them, since $\beta$ uses the
  printed counts. Table XXII's combining rule for its three multi-$\sigma$ cells
  is likewise unprinted, and the geometric mean is stated as a reconstruction.
- **Origins cited but not consulted:** Freundlich, *Kapillarchemie* (1909) — read
  only through Langmuir's restatement and Table VII's $q_F$ column; Eucken
  (1914); Bakker (1915); Haber (1914); Knudsen (1910); Meyer's *Kinetic Theory of
  Gases* (1899); and Langmuir's own Parts I and II. Nothing here derives from any
  of them.
- **`J1.3` is read but not loaded.** No dataset of `J1.3` is loaded and no number
  of `J1.3` is retyped; the §6.8 identity is verified against Langmuir's own
  eq. (30).
- **This page is not a rate.** Adsorption breakthrough and the linear-driving-force
  law are `J1.5`; a Langmuir equilibrium drops into that page's isotherm slot.

## One thing this page checked for somebody else

Rawlings & Ekerdt, *Chemical Reactor Analysis and Design Fundamentals* (2nd edn),
book page 443, state that the Danckwerts boundary conditions "were derived at
least 45 years prior to Danckwerts in a classic paper by Langmuir [22]".
**That paper is not this one.** Their reference [22] is *Langmuir, I., "The
velocity of reactions in gases moving through heated vessels and the effect of
convection and diffusion", J. Am. Chem. Soc.* **30**(11), 1742–1754 (1908) — read
from their own bibliography, and 1953 − 1908 = 45 exactly. The 1918 paper on disk
contains no flow, no convection and no boundary condition of any kind: the words
"convection" and "flow" do not occur in its forty-three pages, "diffusion" occurs
three times and never in a transport equation, and every apparatus in it is a
sealed static bulb. The claim bears on `A2.1`/`A2.2` and **cannot be settled from
this file**; settling it needs JACS **30**(11) 1742, which is not on disk.
Recorded on the case yaml as an acquisition target. `A2.1` and `A2.2` were not
touched.
