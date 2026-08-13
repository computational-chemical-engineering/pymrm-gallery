# B3.3 — Bhatia & Perlmutter's random pore model, Part I (isothermal, kinetic control)

From Bhatia, S. K. & Perlmutter, D. D., *A random pore model for fluid–solid
reactions: I. Isothermal, kinetic control*, AIChE Journal **26**(3), 379–386
(May 1980), [doi:10.1002/aic.690260308](https://doi.org/10.1002/aic.690260308).
Read from the 8-page Wiley scan at its native 300 ppi (`pdfimages -list`:
CCITT-G4 bilevel, 300 × 300 ppi on every page); every equation and every
constant was transcribed from cropped native renders, and the text layer was
used only as a search index.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page; runs in about 8 s |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata |
| `data/bhatia-1980-printed-scalars.csv` | the 47 numbers the paper prints outside its equations |
| `data/random-pore-reference.csv` | the page's own model evaluations, for regression |

## The scope decision, and why it is the first result

**The paper contains no table.** All eight pages were rendered at 300 ppi and
read column-block by column-block: eight **figures** (Figs. 1–8), a NOTATION
list, a LITERATURE CITED list, nothing else tabular. Second, independent
search: the text layer grepped case-insensitively for `table` returns two
lines, both inside words (`suitable`, p. 385; a scan-garbled `surfaces`,
p. 380) — never a caption, never a cross-reference.

The paper's **only** comparison with experiment is **Figure 8**, the char data
of Hashimoto et al. (1979), which is not on disk here. It is **not digitised**:
no point, curve or axis coordinate is taken from it or from any of the other
seven figures. **This page therefore does not establish that the random pore
model describes any real solid**, and says so in its title cell, in §4 and in
Reuse.

What separates this case from the figure-only ones is that the paper prints
its Figure 8 fit results **in prose** on p. 385 — $S_o = 520$ m²/cm³ for both
chars, $\psi = 6.9$ and $\psi = 13.7$ — and those three numbers are transcribed
from that sentence. Everything derived from them is a *prediction* of the model
at those parameters, never a validation of it.

## What the page establishes

- **The derivation chain, symbolically.** Fifteen `sympy` identities, all
  exactly zero, from eqs. (17)–(30) to (31)–(34), (36) and (37). Two are worth
  naming: eq. (37)'s radical is the **perfect square** $(1+\psi\tau/2)^2$, so
  the paper's $S/S_o$ **is** $dX/d\tau$; and Petersen's eqs. (3)–(4) put his
  surface maximum at $\epsilon = 1/2$ for **every** $L$ and $r_p$, which is the
  whole content of eq. (34).
- **The printed structural claims, as numbers.** From the paper's own
  squared-deviation integral $I$: best-fit $m$ = **1.0000000000** at $\psi = 0$
  and **0.4905224** at $\psi = 2$ (printed: $0.49 \le m \le 1$), and
  **0.6638850** at $\psi = 1$, 0.417 % below two-thirds (printed: "very close
  to two-thirds"). $X_M \to 1-e^{-1/2} = 0.3934693$ (printed: $< 0.393$).
  $\tau_{1/2}$ at $\sigma = 100$ within **2.9 %** of its plateau, and spreading
  only **7.2 %** across $0 \le \psi \le 100$ at $\sigma = 0.25$ — the paper's
  two prose claims about Figure 2.
- **Both special-case reductions the abstract claims.** Bhatia's eq. (40) *is*
  the published `B3.1`'s Yagi–Kunii eq. (6) in its reaction-control limit, to
  **1.1e-16** — two transcriptions of two different papers, three months apart.
  And the approach to Petersen (1957) is measured, not asserted: the two
  conversion–time curves stay within **0.0218** in $X$ up to $X = 0.75$ and
  first part by two conversion points at **$X = 0.7423$**, which is the paper's
  "good agreement up to about 75 % conversion" to three digits — **under the
  convention p. 384 prints**, quoted on the page: *"If $\epsilon_o = 0.3$ and
  $L_o = 3.14 \times 10^6$ are chosen, the results are the same as the prior
  curve B."* The Fig. 7 caption prints three numbers for a two-parameter
  Petersen structure, so a convention has to be chosen; this one is the
  paper's, and the whole band is printed beside it rather than one
  representative of it — the other two pairings, $(\epsilon_o, S_o)$ and
  $(S_o, L_o)$, give 0.7599 and 0.7687, so the dependence is at the *second*
  digit and the reader can see which digits belong to the paper's choice. The
  two curves also cross **twice** — Petersen
  starts 0.66 % faster, so curve A trails until $X = 0.0232$. The surfaces
  part much earlier (−21.7 % at $X = 0.75$ on the Petersen-consistent Fig. 4
  set, −21.0 % on curves A and B's own), and the reason is computed: the
  two models use **different overlap laws**, $1-e^{-V_E}$ against
  $3u^2-2u^3$, differing at order $V_E^{3/2}$, so no parameter choice makes
  them identical.
- **One thing the paper's own numbers say that it does not.** Figure 7's curve
  C is called "no longer consistent with the requirements of Petersen's model".
  It is worse than that: eq. (3) has a maximum in $r_p$, so no Petersen
  structure of length $3.14\times10^6$ cm/cm³ can have a surface above
  **2720 cm²/cm³**, and curve C asks for 12 500 — a factor **4.60** beyond it,
  with **no real root at all**. The same arithmetic reconstructs the curve-D
  length the paper says it computed internally and never prints,
  $7.747\times10^7$ cm/cm³.
- **The root of eq. (7), selected by construction.** The paper writes "the
  solution to the cubic equation"; it has three real roots. §2.1 shows eq. (7)
  is exactly $G = 3/(2u_o)$, which picks one — and eq. (6) then reaches
  $X = 1$ exactly at $s = 2G/3-1$, while the middle root puts complete
  gasification at a negative time — and injecting it moves §6.3's conversion
  comparison from 0.0218 to 36.4, a factor **1672**: the 5th largest absolute
  shift among the 53 break rows, though only the 22nd largest relative one, and
  the page prints both ranks rather than calling it the largest. It does **not**
  move `petersen_X_endpoint_dev`, for the reason given below.

## Three printed defects, reported and not repaired

1. **"This match at $m = 1$"** [*sic*, p. 383] — in a sentence whose
   predecessor says the grain model matches the $\psi = 1$ results, whose
   dashed line is labelled $m = 2/3$, and whose "above-mentioned correspondence
   as $\psi \to 0$" is *already* the $m = 1$ case. The least squares settles
   it: at $\psi = 1$ the best $m$ is 0.6639 and $m = 1$ misfits **5.34×**
   worse. Which repair was intended cannot be decided from the document, so
   neither is asserted.
2. **The optimal-porosity worked example.** The printed $\epsilon_o = 0.1$ is
   not a root of the paper's own eq. (35); the root is **0.0852877**, residual
   −4.2e-3. A defect of precision, not substance — the optimum is flat enough
   that 0.1 costs 0.0028 % of the maximum surface — and reported anyway,
   because eq. (35) is an equation. The example's companion threshold does
   check out: $1+\ln 0.7 = 0.6433$ against the printed "≤ 0.64 violates (36)".
3. **Park and Levenspiel (1976)** in the p. 380 text against **(1975)** in the
   LITERATURE CITED, the paper's only Park reference.

Recorded without further comment: "not clearly indentified" [*sic*, p. 382],
"This research was support by the U.S. Department of Energy" [*sic*, p. 385],
and "stoichimetric" [*sic*] in the NOTATION list two lines below
"stoichiometric".

## The checks that can fail, and the ones that cannot

Every metric has a break-table row that moves it (53 rows, 45 metrics), and the
coverage map is assembled from the table itself and asserted key-for-key
against `agreement.json` in both directions. The six metrics below
`check_agreement.py`'s `ABS_FLOOR = 1e-12` are identified **from the measured
values**, not typed, and each has a named above-floor companion whose quoted
numbers are interpolated from the row they refer to.

**"It moved" is not enough.** The coverage assert requires the move to clear a
stated noise floor (1e-9), not merely to be non-zero — a row that shifts a
2e-16 metric to 2e-15 has replaced float noise with float noise. One row here
was exactly that, and the page prints it rather than quietly replacing it: the
wrong root of eq. (7) **cannot** break `petersen_X_endpoint_dev`, because
$X(s = 2G/3-1) = 1$ is an exact identity for *every* root of eq. (7) —
$(4/27)\epsilon G^3 = G-1$ is precisely what collapses the bracket to
$1/\epsilon - 1$. The metric is kept and its real content named: it establishes
that eqs. (6) and (7) are mutually consistent *as transcribed*, and says
nothing about which root was selected. Its row now mis-transcribes eq. (7)
itself, $4/27$ as $2/27$, which moves it to $1/(1-\epsilon_o) = 1.4286$ for any
root of the misread cubic. The smallest move in the table and the weakest row
in relative terms are printed too — `nrms_at_psi_m_star_zero` moves only 3.4 %
under its injected defect, below `check_agreement.py`'s `REL_TOL = 5 %`, so
that one defect would pass the regression suite.

**Three candidate break rows were found not to work, and are reported with a
measured witness each rather than dropped.** Flipping the sign inside
eq. (31)'s radical gives a function real only for $w > e^{-1/\psi}$ and
strictly increasing there, so there is no maximum to root-find at all.
Inverting eq. (35)'s log argument gives residuals of one sign at both ends of
$(0,\epsilon)$, so no root. Halving eq. (20)'s exponent puts Avrami below
Petersen for every $u > 0$, so there is no crossing.

**One metric is structural**, `sym_zero_identity_count`: it counts `sympy`
identities, so it verifies the paper's algebra rather than this page's
numerics, no runtime defect can move it, and it is cited as evidence for
nothing numerical.

**Two identities are narrower than they look, and the page says so.** §7.2's
Petersen routes (c) and (d) share no code but *do* share algebra — substituting
$G = 3/(2u_o)$ turns route 1 into route 2 — so what they test is the selection
of the root of eq. (7) and the transcription of eq. (4), not eq. (6) itself.
And the cross-page identity with `B3.1` is taken at $\omega = 0$, where every
$\omega$ term in eq. (6) drops out **and $\gamma$ cancels identically** (the
numerator is $\gamma x$, the denominator $\gamma$) — so it can see neither the
three $\omega$ coefficients nor the value of $\gamma$, and the $\gamma$ range
printed beside it cannot fail and is labelled as such. It also cannot see a
*consistent* slip in the shrinking-core exponent, which is written once and
used on both sides. What it does establish is that eq. (40) and
eq. (6)$|_{\omega\to0}$ are the same normalised law, and the two axes it can
see get a break row each: the exponent on $x$, and the ratio of the two
$\gamma$ coefficients. `B3.1`'s published `eq6_to_reaction_limit = 3.848988e-10`
is reproduced here to 3.5e-06 at that page's own $\omega = 10^{-9}$, which
shows the 3.85e-10 is entirely its regularisation and not a transcription
difference.

## Traps for anyone reusing this

- $\psi$ carries $(1-\epsilon_o)$; eq. (32) carries $\psi\tau/4$, not
  $\psi\tau/2$. Both have break rows.
- Eq. (3) is **not monotone** in $r_p$ — two roots for any attainable $S$, none
  at all above $S_{\max}(L)$.
- Eq. (32) cannot be inverted naively at large $\psi$: the exponential
  underflows and `newton` gets an exactly singular Jacobian. The page solves
  the log form, with $\tau$ through a sigmoid so it cannot leave $(0,\sigma)$.

## Scope, honestly

Kinetic control throughout — the paper's own title. **Part II** (non-isothermal
operation and diffusion control) is not on disk, was not consulted, and nothing
here depends on it. Petersen (1957) and Avrami (1940) are read **only** as this
paper reprints them, so every statement about "Petersen's model" here is a
statement about that reprint. `B3.2` is the neighbouring page where intergrain
diffusion is present and matters; `B3.1` is the shrinking core this model
collapses onto.
