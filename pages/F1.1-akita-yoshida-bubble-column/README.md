# F1.1: Akita and Yoshida, gas holdup and kLa in bubble columns

Gas holdup and volumetric mass-transfer coefficient correlations from Akita &
Yoshida, *Ind. Eng. Chem. Process Des. Dev.* **12**(1) 76-80 (1973),
[doi:10.1021/i260045a015](https://doi.org/10.1021/i260045a015), and the bubble
size, interfacial area and liquid-side coefficient of their sequel, **13**(1)
84-91 (1974), [doi:10.1021/i260049a016](https://doi.org/10.1021/i260049a016).

The page checks the two papers against their own tables, one digitised figure and
each other: the 1974 kL times the 1974 area against the 1973 kLa at the same
holdup.

## Contents

| file | what it is |
| --- | --- |
| `index.py` | the page, in jupytext percent format. Edit this, not the notebook |
| `index.ipynb` | the paired, executed notebook; `jupytext --sync` keeps the two in step |
| `agreement.json` | every number the page states, written by the Claims cell and checked by CI |
| `data/akita-yoshida-1974-table2-properties.csv` | liquid properties, 1974 Table II |
| `data/akita-yoshida-1974-table4-bubble-sizes.csv` | the counted bubble size distribution of one run, 1974 Table IV |
| `data/akita-yoshida-1973-table2-ranges.csv` | ranges and dimensionless groups of the mass-transfer runs, 1973 Table II |
| `data/akita-yoshida-1974-fig12.csv` | calculated against observed kLa, digitised from 1974 Figure 12 |
| `data/akita-yoshida-printed-numbers.csv` | every published constant and printed result the page uses |
| `meta.yaml` | page metadata, validated against `models.yaml` |

Each CSV has a `.meta.yaml` provenance sidecar beside it.

## Things to know before reusing this

**This is a correlation, not a discretised model.** The functions `holdup_eq11`,
`kla_eq18`, `dvs_eq19`, `kl_eq32` and `a_eq24` in `index.py` take SI arguments
and the dictionary of published constants; no pymrm operator appears. pymrm
supplies `pymrm.checks.find_roots` for the implicit holdup.

**The two papers' kLa are not one correlation.** The 1974 product of kL and area
is 1.071 <!-- claim: kla_product_geo_ratio --> times the 1973 kLa at the same
holdup, as a geometric mean over the matched 1973 Table II points.

**Eq. 23 is a power-law stand-in for the holdup correlation** that the 1974 paper
uses up to the holdup limit it states. At that limit eq. 11 needs
1.280 <!-- claim: eq23_ratio_at_eps_max --> times the velocity eq. 23 gives for the
same holdup. Use eq. 11.

**One Table IV entry disagrees with its counts.** The third-moment entry of the eleventh size class is
printed 25.54 <!-- claim: t4_worst_dmu3_printed --> where the counts give 23.55 <!-- claim: t4_worst_dmu3_computed -->;
the printed column total agrees with the computed value, not with
the printed entry.
