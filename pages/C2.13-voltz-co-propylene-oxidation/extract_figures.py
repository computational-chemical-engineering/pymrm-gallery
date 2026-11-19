#!/usr/bin/env python3
"""Regenerate data/voltz-1973-fig19-fig20.csv from the source PDF.

The PDF is not in the repository and must not be (AGENTS.md). Point this at your
own lawfully obtained copy:

    python extract_figures.py <your copy>.pdf

Requires poppler's ``pdftoppm`` on PATH, plus numpy, scipy and pillow.

Method
------
Only the experimental markers are extracted; the solid lines are the authors'
own calculations and are left alone. Each marker was located by eye on a
600 dpi render of the figure (the approximate positions below, in the pixel
coordinates of the crop), and assigned to its curve A to G by its glyph
(A inverted triangle, B circle, C cross, D triangle, E diamond, F circle,
G square; B and F are both circles and lie two decades apart in space
velocity). The position is then refined automatically: within +/- 10 px of the
eye reading, the centre is the location of the largest ink count in a 23 px
box, which is the glyph rather than a thin calculated line passing near it.
Where two glyphs touch (eye readings closer than 30 px, two pairs on each
figure), the box would pull both onto the darker one, so those keep the eye reading.

The authors' calculated lines are read at one conversion level each (50 % on
Figure 19, 40 % on Figure 20): along the pixel row of that level, every run of
ink is listed, and the run that belongs to each line is named in LINE_RUNS
after inspection (run 0 is the y-axis tick; on Figure 20 the 40 % row is used
because at 50 % the C and D lines pass through markers, and D is left out at
40 % for the same reason). These positions are used only to test the
reconstruction of the authors' calculation, never as data.

Axes are calibrated on the printed tick marks of the crop, located by column
and row sums of the ink: the logarithmic x axis from the ticks at 0.1, 0.2,
1.0, 2.0 and 20, the y axis from the nine ticks at 10 to 90 %.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import uniform_filter

DPI = 600
# (pdf page, crop box (col_lo, row_lo, col_hi, row_hi)) in the 600 dpi render
CROP = {19: (6, (2650, 1220, 4900, 2360)), 20: (7, (450, 320, 2560, 1420))}
# Axis calibration in crop pixels: column of SV = 0.1, pixels per decade,
# row of 0 % and pixels per 1 %.
CAL = {19: (249.0, 585.3, 950.3, 8.78), 20: (188.0, 593.0, 960.5, 8.95)}

# Eye readings in crop pixels (x, y), per curve.
EYE = {
    19: {
        "A": [(290, 81), (299, 161), (310, 219), (346, 273), (405, 431),
              (475, 578), (614, 729), (713, 794)],
        "B": [(492, 132), (512, 214), (539, 398), (608, 603), (741, 768),
              (835, 821)],
        "C": [(648, 82), (676, 205), (723, 347), (794, 484), (845, 574),
              (894, 664), (1031, 790), (1167, 862)],
        "D": [(785, 150), (852, 258), (898, 362), (963, 461), (1038, 593),
              (1091, 655), (1196, 767), (1334, 853)],
        "E": [(938, 116), (997, 246), (1104, 435), (1177, 600), (1202, 656),
              (1344, 779)],
        "F": [(1042, 168), (1074, 245), (1079, 257), (1114, 335), (1170, 450),
              (1175, 475), (1190, 538), (1258, 611), (1297, 681), (1335, 695)],
        "G": [(1141, 183), (1166, 232), (1194, 281), (1256, 423), (1291, 486),
              (1298, 574), (1349, 592), (1398, 660), (1476, 740)],
    },
    20: {
        "A": [(230, 221), (241, 456), (250, 548), (285, 654), (333, 810),
              (413, 907), (551, 942), (660, 947)],
        "B": [(439, 233), (456, 376), (483, 573), (549, 744), (684, 881),
              (783, 907)],
        "C": [(620, 322), (667, 503), (737, 635), (779, 710), (836, 789),
              (968, 863), (1102, 923)],
        "D": [(734, 184), (802, 361), (913, 603), (987, 691), (990, 710),
              (1050, 789), (1148, 830), (1293, 899)],
        "E": [(847, 503), (886, 262), (942, 451), (1038, 648), (1128, 755), (1153, 793),
              (1293, 860)],
        "F": [(995, 253), (1036, 398), (1075, 471), (1127, 581), (1127, 610),
              (1210, 721), (1257, 775), (1292, 785)],
        "G": [(1097, 264), (1120, 350), (1148, 414), (1214, 554), (1250, 612),
              (1303, 707), (1358, 762), (1430, 818)],
    },
}

# Calculated lines: (conversion level %, {curve: index of the ink run in that row})
LINE_RUNS = {19: (50, dict(zip("ABCDEFG", range(1, 8)))),
             20: (40, {"A": 1, "B": 2, "C": 3, "E": 7, "F": 8, "G": 10})}

WIN = 23      # glyph box, px
SEARCH = 10   # search radius around the eye reading, px


def render(pdf: Path, page: int, tmp: Path) -> np.ndarray:
    stem = tmp / f"p{page}"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(page), "-l", str(page),
                    "-gray", "-png", str(pdf), str(stem)], check=True)
    png = next(tmp.glob(f"p{page}*.png"))
    return np.asarray(Image.open(png).convert("L"))


def refine(ink: np.ndarray, x: int, y: int) -> tuple[float, float]:
    dens = uniform_filter(ink.astype(float), WIN)
    win = dens[y - SEARCH:y + SEARCH + 1, x - SEARCH:x + SEARCH + 1]
    j, i = np.unravel_index(np.argmax(win), win.shape)  # densest glyph box
    return x - SEARCH + i, y - SEARCH + j


def row_runs(ink: np.ndarray, row: int, col0: int, col1: int) -> list[float]:
    """Centre columns of the ink runs along a 3 px band about ``row``."""
    seg = ink[row - 1:row + 2, col0:col1].any(axis=0)
    edges = np.flatnonzero(np.diff(np.r_[0, seg.astype(int), 0]))
    return [col0 + 0.5 * (a + b - 1) for a, b in zip(edges[::2], edges[1::2])]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", type=Path)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent / "data" / "voltz-1973-fig19-fig20.csv")
    args = ap.parse_args()
    rows, lines = [], []
    with tempfile.TemporaryDirectory() as t:
        for fig, (page, box) in CROP.items():
            img = render(args.pdf, page, Path(t))
            crop = img[box[1]:box[3], box[0]:box[2]]
            ink = crop < 128
            x0, dec, y0, pct = CAL[fig]
            every = np.array([p for pts in EYE[fig].values() for p in pts], float)
            for curve, pts in EYE[fig].items():
                for (x, y) in pts:
                    dist = np.hypot(every[:, 0] - x, every[:, 1] - y)
                    touching = np.sort(dist)[1] < 30
                    xr, yr = (x, y) if touching else refine(ink, x, y)
                    sv = 10 ** (-1 + (xr - x0) / dec)
                    conv = (y0 - yr) / pct
                    rows.append((fig, curve, sv, conv, xr - x, yr - y))
            level, which = LINE_RUNS[fig]
            runs = row_runs(ink, int(round(y0 - pct * level)), int(x0) + 8,
                            crop.shape[1] - 600)
            for curve, i in which.items():
                lines.append((fig, curve, level, 10 ** (-1 + (runs[i] - x0) / dec)))
    species = {19: "CO", 20: "C3H6"}
    with open(args.out, "w") as f:
        f.write("# Voltz, Morgan, Liederman & Jacob (1973), Ind. Eng. Chem. Prod. Res.\n"
                "# Dev. 12, 294-301, Figures 19 (CO) and 20 (C3H6), 550 F, synthetic gas.\n"
                "# Experimental markers only; regenerate with extract_figures.py.\n")
        f.write("figure,species,curve,space_velocity,conversion_pct\n")
        for fig, curve, sv, conv, dx, dy in rows:
            f.write(f"{fig},{species[fig]},{curve},{sv:.4g},{conv:.1f}\n")
    with open(args.out.parent / "voltz-1973-calculated-lines.csv", "w") as f:
        f.write("# Voltz et al. (1973), Figures 19 and 20: space velocity at which the\n"
                "# authors' calculated (solid) line crosses one conversion level.\n"
                "# Regenerate with extract_figures.py.\n")
        f.write("figure,species,curve,conversion_pct,space_velocity\n")
        for fig, curve, level, sv in lines:
            f.write(f"{fig},{species[fig]},{curve},{level},{sv:.4g}\n")
    shifts = np.array([[r[4], r[5]] for r in rows])
    print(f"{len(rows)} markers; refinement shift max {np.abs(shifts).max():.0f} px")


if __name__ == "__main__":
    main()
