#!/usr/bin/env python3
"""Stage the Figure 9 marker extraction for Kunii & Levenspiel (1968), page E2.1.

Run locally, against a lawfully obtained copy of the paper. Nothing this script
reads or writes may enter the repository: `pdftoppm` output and the overlay are
page images, and the gallery does not commit those.

    python extract_figure9.py ~/papers/pymrm-gallery/i260028a001.pdf /tmp/e21

Writes to the output directory:
    fig9_candidates.csv   detected centres in data coordinates, with shape scores
    fig9_overlay.png      candidates drawn on the 600 dpi render, for visual review
    fig9_contact.png      44 px patches around every candidate, for auditing

WHAT THE MAINTAINER IS BEING ASKED
----------------------------------
Figure 9 is conversion (1 - X_A, log axis, 1.0 down to 0.04) against the
dimensionless rate group K_m (linear, 0 to 8), for Kobayashi & Arai's ozone bed.
Three series are distinguished only by glyph:

    (dotted circle)  u_0 = 6.6 cm/s    u_0/u_mf = 3.14
    (filled circle)  u_0 = 9.9 cm/s    u_0/u_mf = 4.71
    (open circle)    u_0 = 13.2, 16.5, 20 cm/s   u_0/u_mf = 6.28, 9.55

and the paper's claim under test is that the fitted bubble size rises with
velocity (d_b = 3.7, 4.2, 5.0 cm for the three series respectively).

The extraction splits cleanly in two:

  * K_m greater than about 1.8 - markers are isolated, the ring template locks on,
    and the three glyphs separate on the interior ink fraction (filled: hole ink
    ~1.0; dotted: hole ~0.1-0.2 with a solid 5 px centre; open: hole < 0.1 and no
    centre). Roughly 14 markers. This region is where the three model curves are
    far apart, so it carries essentially all of the discriminating power.

  * K_m below about 1.8 - roughly 20 markers overlap into a connected black mass
    in which individual centres cannot be located at all, let alone glyphs. See
    fig9_contact.png patches 0-15. This is not a classifier problem; the ink is
    merged. Those markers are unresolvable and would have to be reported as such.

QUESTION FOR THE MAINTAINER: is the right-hand subset alone worth publishing as
the experimental comparison for E2.1, given that (a) d_b was fitted per series by
the authors, so the comparison tests the fit rather than a prediction, and (b)
about 20 of the roughly 34 markers would be dropped as unresolvable? If yes, the
page gains an experimental tier; if no, E2.1 stays tier 6 on the three worked
appendices, which is how it is currently written.

AXIS CALIBRATION (established, and independent of the marker question)
----------------------------------------------------------------------
On the 600 dpi render of PDF page 9, with the sub-image taken as
rows 200:1700, columns 2550:4900, the plot frame is

    top row 234.5, bottom row 1290.5, left column 766.5, right column 1739.5

A least-squares fit to the printed ticks gives

    row   = 236.660 - 750.532 * log10(1 - X_A)     residuals <= 5 px over 15 ticks
    col   = 764.544 + 121.475 * K_m                residuals <= 2 px over 9 ticks

and the frame corners recover 1.007 and 0.0394 against the nominal 1.0 and 0.04,
so the axis limits are exactly 1.0 and 0.04 as labelled.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

# sub-image window and frame, in pixels of the 600 dpi render of PDF page 9
WIN = (200, 1700, 2550, 4900)          # row0, row1, col0, col1
FRAME = dict(top=234.5, bottom=1290.5, left=766.5, right=1739.5)
Y_FIT = (236.660, -750.532)            # row = a + b*log10(1-X_A)
X_FIT = (764.544, 121.475)             # col = a + b*K_m


def render(pdf: Path, out: Path) -> np.ndarray:
    subprocess.run(["pdftoppm", "-r", "600", "-f", "9", "-l", "9", "-png",
                    str(pdf), str(out / "p")], check=True)
    im = Image.open(next(out.glob("p-*.png"))).convert("L")
    a = np.asarray(im)
    return a[WIN[0]:WIN[1], WIN[2]:WIN[3]]


def ann(r0, r1, rmax=22):
    y, x = np.mgrid[-rmax:rmax + 1, -rmax:rmax + 1]
    d = np.hypot(x, y)
    return ((d >= r0) & (d <= r1)).astype(float)


def envelope(Km, db=5.0):
    """The topmost printed model curve, 1 - X_A for d_b = 5.0 cm.

    No experimental marker can sit above it by more than the line width, so this
    is a principled filter for the glyphs of the annotation block ("Calcd. by
    Eqs. 54 and 56", "d_b = 5.0 cm", "cm"), which the ring template also matches.
    Equations 2, 3, 10, 12, 15, 17, 45, 49, 50 and 51 for appendix C's bed.
    """
    g, u0, umf, emf, em, Lm, D, alpha = 980.0, 13.2, 2.1, 0.50, 0.45, 34.0, 0.204, 0.47
    ubr = 0.711 * np.sqrt(g * db)
    ub = u0 - umf + ubr
    delta = (u0 - umf) / ub
    Kbc = 4.5 * (umf / db) + 5.85 * (D ** 0.5 * g ** 0.25 / db ** 1.25)
    Kce = 6.78 * np.sqrt(emf * D * ub / db ** 3)
    gc = (1 - emf) * (3.0 * (umf / emf) / (ubr - umf / emf) + alpha)
    ge = (1 - emf) * (1 - delta) / delta - gc
    Kr = Km * u0 / ((1 - em) * Lm)
    brack = 1.0 / (Kr / Kbc + 1.0 / (gc + 1.0 / (Kr / Kce + 1.0 / ge)))
    return np.exp(-(1.0 / (1 - emf)) * (u0 / ubr) * Km * brack)


def detect(sub: np.ndarray):
    ink = (sub < 128).astype(float)
    mf = lambda m: ndimage.correlate(ink, m / m.sum(), mode="constant")
    # marker geometry measured on an isolated glyph: outer ring r = 11-15.5,
    # white hole r = 4-9.5, centre dot r <= 2.5
    ring, hole, dot = mf(ann(11.0, 15.5)), mf(ann(4.0, 9.5)), mf(ann(0, 2.5))
    cand = ring > 0.90
    T, B, L, R = (int(FRAME["top"]), int(FRAME["bottom"]),
                  int(FRAME["left"]), int(FRAME["right"]))
    cand[:T + 12] = cand[B - 12:] = False
    cand[:, :L + 12] = cand[:, R - 12:] = False
    lab, n = ndimage.label(cand)
    out = []
    for i in range(1, n + 1):
        ys, xs = np.where(lab == i)
        cy, cx = float(ys.mean()), float(xs.mean())
        iy, ix = int(round(cy)), int(round(cx))
        h, d = float(hole[iy, ix]), float(dot[iy, ix])
        shape = "filled" if h > 0.8 else ("dotted" if d > 0.8 else "open")
        Km = (cx - X_FIT[0]) / X_FIT[1]
        one_minus_X = 10.0 ** ((cy - Y_FIT[0]) / Y_FIT[1])
        out.append(dict(row=cy, col=cx, Km=Km, one_minus_XA=one_minus_X,
                        hole=h, dot=d, shape=shape, npix=len(ys),
                        resolved="yes" if Km > 1.8 else "no",
                        suspect_text="yes" if one_minus_X > 1.05 * envelope(Km) else "no"))
    out.sort(key=lambda r: r["col"])
    return out, ink


def main():
    pdf = Path(sys.argv[1]).expanduser()
    out = Path(sys.argv[2]).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    sub = render(pdf, out)
    cands, _ = detect(sub)

    cols = ["Km", "one_minus_XA", "shape", "resolved", "suspect_text", "hole", "dot", "row", "col", "npix"]
    lines = [",".join(cols)]
    for c in cands:
        lines.append(",".join(f"{c[k]:.5g}" if isinstance(c[k], float) else str(c[k])
                              for k in cols))
    (out / "fig9_candidates.csv").write_text("\n".join(lines) + "\n")

    rgb = Image.fromarray(sub).convert("RGB")
    dr = ImageDraw.Draw(rgb)
    colour = {"filled": (220, 30, 30), "dotted": (20, 120, 220), "open": (20, 160, 40)}
    for k, c in enumerate(cands):
        x, y = c["col"], c["row"]
        dr.ellipse([x - 20, y - 20, x + 20, y + 20], outline=colour[c["shape"]], width=3)
        dr.text((x + 22, y - 22), str(k), fill=colour[c["shape"]])
    rgb.save(out / "fig9_overlay.png")

    P, ncol = 44, 8
    nrow = (len(cands) + ncol - 1) // ncol
    sheet = Image.new("L", (ncol * (2 * P + 8), nrow * (2 * P + 22)), 255)
    ds = ImageDraw.Draw(sheet)
    for k, c in enumerate(cands):
        iy, ix = int(round(c["row"])), int(round(c["col"]))
        r, col = divmod(k, ncol)
        sheet.paste(Image.fromarray(sub[iy - P:iy + P, ix - P:ix + P]),
                    (col * (2 * P + 8) + 4, r * (2 * P + 22) + 18))
        ds.text((col * (2 * P + 8) + 4, r * (2 * P + 22) + 2),
                f"{k} {c['shape'][:4]} h={c['hole']:.2f}", fill=0)
    sheet.save(out / "fig9_contact.png")

    nres = sum(c["resolved"] == "yes" for c in cands)
    print(f"{len(cands)} candidates, {nres} in the resolvable region K_m > 1.8")
    print(f"wrote {out}/fig9_candidates.csv, fig9_overlay.png, fig9_contact.png")
    print("NOTE: the overlay and contact sheet are page images. Do not commit them.")


if __name__ == "__main__":
    main()
