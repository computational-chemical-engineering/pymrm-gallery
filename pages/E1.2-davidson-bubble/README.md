# E1.2 — The Davidson bubble

A single bubble in an incipiently fluidised bed is a Laplace problem with a
uniform-pressure hole in it. Solve it once and you get the cloud, the
through-flow, and — because both come from the same field — the reason two
constants that look unrelated in the reactor literature are not.

- **Structures:** `S3` (1D steady BVP)
- **Runtime:** ~9 s
- **Tier:** 6. Nothing here is compared with a measurement.

## Sources, and the distinction that matters

**Kunii, D. & Levenspiel, O.**, *Bubbling bed model for kinetic processes in
fluidized beds*, Ind. Eng. Chem. Process Des. Dev. **7**(4) 481–492 (1968),
[doi:10.1021/i260028a001](https://doi.org/10.1021/i260028a001) — **the only text
read.** It prints the results as its equations 3, 9, 10, 11 and 17 and, what makes
the reprint carry the case, *evaluates* each of them with the arithmetic shown and
the answer printed in three worked appendices.

**The attribution, at its actual strength.** The name Davidson appears in the body
of that paper exactly **once**: on journal page 490, inside a parenthesis listing
the model's simplifying assumptions — *"many simplifying assumptions were made in
developing this model (the Davidson bubble, bed with single size of bubbles,
…)"*. That is the whole of it. There is **no** Davidson & Harrison entry in their
literature cited (journal page 492), where `Davidson, J. F.` appears only as a
co-author of Orcutt, Davidson & Pigford (1958); equations 3 and 9 and the `0.711`
carry no citation at all; and every derivation is deferred to a companion paper
(Kunii & Levenspiel, *I&EC Fundamentals* **7** 446, 1968b) that is not on disk.

The monograph usually named as the origin of these results was **not consulted**
and is the source of nothing here — and even the knowledge that it is the usual
citation is general background rather than something traced to a document read for
this page, so its bibliographic detail is deliberately absent from the metadata.
**The derivation on the page is done here from scratch.**

The slug is `davidson-bubble`, not `davidson-harrison-bubble`: nothing read here
attributes the model to two authors.

## Contents

| file | what it is |
| --- | --- |
| `index.ipynb` | the page |
| `build_page.py` | regenerates `index.ipynb`. Edit this, not the notebook |
| `meta.yaml` | page metadata |

**No `data/` directory.** The printed appendix numbers this page checks against
are already transcribed, sidecar and all, in
`pages/E2.1-kunii-levenspiel-bubbling-bed/data/`, and are loaded **cross-page**
with `load_data(..., page="E2.1-kunii-levenspiel-bubbling-bed")`. One
transcription in the repository, no chance of two copies drifting. Every row used
was independently re-read here on a 600 dpi render before use. **This page must
be integrated while `E2.1` remains published** — the same dependency `F1.3` has
on `F1.4` and `A4.2` has on `A4.9`.

## The model in one paragraph

In the frame of the bubble, the solids stream down past it at `u_br` in inviscid
potential flow, and the interstitial gas percolates through those solids at
`u_f = u_mf/eps_mf` by Darcy's law. The bubble contains gas of negligible
viscosity and weight, so its surface is at **uniform pressure**. Both the solids
potential and the gas potential are then harmonic with the same `cos(theta)`
angular dependence, so both reduce to the **same 1D radial ODE**

```
(1/r^nu) d/dr( r^nu dh/dr ) - nu h / r^2 = 0
```

with `nu = 2` for a spherical bubble and `nu = 1` for the 2-D one. One operator,
two boundary conditions at `r = R`: `h = 0` for the gas (uniform pressure),
`dh/dr = 0` for the solids (nothing crosses).

## Validation, in the order the brief ranks it

**Route 1 — worked examples with printed intermediates.** All three of Kunii and
Levenspiel's appendices evaluate one of these quantities with the arithmetic
shown. Appendix C's `gamma_c` line is a fully substituted evaluation of
**equation 9** with a printed answer, which makes the cloud volume a worked
example rather than a formula taken on trust. **But all eight route-1 values are
already published on `E2.1`, to every digit, from the one shared transcription
this page loads cross-page.** They are the same reading checked twice, not
independent corroboration, and the page says so.

**Route 2 — an internal identity, and this is the new evidence.** The through-flow
computed from the numerical field **rebuilds the `4.5`** that opens equations 10
and 11. That constant was never an input; it falls out of `-f'(R) = 3 u_f` and a
hemisphere integral. The `1.1e-6` is discretisation error, not the strength of the
agreement — the content of the check is binary, and check 9 measures how loudly it
fails.

**No figure was digitised, and none is needed.**

| check | result |
| --- | --- |
| through-flow coefficient vs the printed 4.5 | 1.1e-6 |
| smallest shift of that coefficient under 7 injected operator defects | 33 % |
| cloud radius vs the closed form, n = 1600 | 7.2e-8 |
| observed order of the discretisation | 2.000 |
| cloud radius, root-of-W vs streamline integration | 1.2e-14 |
| 2-D cloud vs its own closed form | 3.9e-8 |
| `K_bc`, `H_bc` rebuilt vs printed 46.5 / 5.44 / 0.0360 (E2.1's numbers) | 0.43 % worst |
| printed appendix values `12.0`, `5.8`, `42.8`, `53.9`, `0.40` (E2.1's numbers) | 1.06 % worst |

**Two reported numbers are weaker than they look, and check 9 says how much.** The
discrete integral identity closes to 2.0e-12, but it is an *algebraic* identity —
the volume-weighted sum of the discrete residuals — and it telescopes for any
consistent operator with any boundary data. It catches a **mismatch** between the
`nu` in `construct_div` and the one in the sink term, and nothing else: it reads
*better* than baseline under a wrong-but-consistent `nu` that makes the
coefficient 77 % wrong, and is unchanged under a mis-signed outer boundary
condition that ruins the field by seven orders of magnitude. `R_c/R` is likewise
blind to a common-mode far-field error. **If you copy this skeleton, watch the
through-flow line, not the conservation line.**

## The finding worth carrying forward

`V_c/V_b = 3/(u_br/u_f - 1)` has a **pole**, not a zero. Below `u_br = u_f` there
is no cloud at all and equation 9 returns a negative number.

**Appendix B of the source paper is on the wrong side of that threshold.** Its
printed `d_b = 0.50` cm, `u_mf = 10` cm/s and `eps_mf = 0.50` give
`u_br/u_f = 0.787` and `V_c/V_b = -14.1`. Nothing in the paper is wrong —
appendix B is a heat-transfer example that needs `gamma_b` and `H_bc` and never
evaluates equation 9 — but the bubble there would have to exceed 0.81 cm for a
cloud to exist, and a reader who carries equation 9 into a bed like that gets
nonsense with no *quantitative* warning: the paper says equation 9 is for "beds
with fast rising bubbles" and never says how fast is fast. The verdict is robust
to the one appendix-B number in doubt — `E2.1` records that its printed
`(1-eps_f)u_b = 8.70` implies `eps_mf = 0.447` rather than the stated `0.50`,
which moves `u_br/u_f` from 0.787 down to 0.704, and nothing below
`eps_mf = 0.635` gives that bed a cloud.

The through-flow has no such threshold: `u_br` never enters the gas potential, so
the `4.5` is valid on both sides of it while the `3` is not. The two constants
are usually quoted together as though they had the same standing.

## The reusable trick

Only two radial modes exist, `h = A r` and `h = B r^-nu`, and the combination

```
dh/dr + (nu/r) h = (nu+1) A
```

is satisfied **exactly by both**. Imposing it as a Robin condition at a finite
outer radius pins the far field without truncating the decaying mode, so a finite
domain carries an infinite one with no truncation error and grid refinement
measures the discretisation alone. That is what makes the observed order of
2.000 mean something.

## What this page is not

It is not `E2.1`. That page takes `u_br`, `V_c/V_b` and `K_bc` as **inputs** and
builds the bubbling-bed reactor model (`S7`) on them; it never solves a flow
field. This page is the `S3` hydrodynamic result those inputs come from.
