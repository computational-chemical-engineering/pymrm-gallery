"""Digitise Figure 8 of Saez and Carbonell (1985), journal p. 60 (PDF page 9).

Usage (from the page directory; the PDF is never stored in the repository):

    python data/digitise.py /path/to/Saez1985-...-AIChEJ31-52.pdf

The page is extracted at its native 300 ppi (``pdfimages -png``, bilevel, no
resampling) into a temporary folder, read, and deleted. Steps:

1. Axes: the frame column and row, and the tick marks at 10^-1 and 10^0 on
   both axes (dark runs of at least 8 px next to the frame); each axis is a
   least-squares line of log10(value) on pixels through three points (the
   frame corner is 10^-2 on both axes).
2. Drawn lines: every dark pixel in the plot area is put on its offset
   D = log10(cal) - log10(exp) from the tick calibration; the solid diagonal
   and the two dashed lines are the centroids of the three peaks of D,
   refined in shrinking windows.
3. Markers: pixels within 0.0085 decades (about 3 px) of a drawn line are
   removed, the rest dilated by 2 px and labelled; each blob larger than 25 px
   is a candidate. Its centre is the centre of its undilated bounding box
   (circles, squares) or the vertex centroid (triangles: the bounding-box centre
   moved by h/6 toward the base). Symbols are assigned by hand from a 3x zoom
   (``HAND`` below); blobs holding several touching markers are ``cluster`` rows,
   with the range of D over their dark pixels, and two touching pairs whose
   members can be told apart are split by hand reads.

It writes ``saez-carbonell-1985-fig8.csv`` and prints the tick residuals, the
line offsets and the marker table the sidecar quotes.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

HERE = Path(__file__).parent
OUT = HERE / "saez-carbonell-1985-fig8.csv"
PDF_PAGE = 9  # journal p. 60

# Plot window (native px) that holds the markers and the three lines.
ROWS, COLS = (230, 950), (255, 985)
LINE_HALF = 0.0085   # decades removed around each drawn line before labelling
MIN_AREA = 25        # px of a marker blob

# Hand assignment, from a 3x zoom: (approximate centre x, y in native px) -> symbol.
# Table 6 (p. 58): open/filled inverted triangle Larkins rings/spheres, open square
# Charpentier et al., open circle Matsuura and Akehata, open triangle Beimesch rings.
HAND = {
    (521, 662): ("filled inverted triangle", "isolated; crossed by the +20% line"),
    (589, 556): ("open inverted triangle", "isolated"),
    (637, 524): ("open square", "isolated"),
    (640, 496): ("filled inverted triangle", "isolated"),
    (676, 573): ("open triangle", "touching an open triangle of the lower cluster"),
    (699, 462): ("filled inverted triangle", "isolated"),
    (728, 420): ("open inverted triangle", "isolated"),
    (749, 499): ("open circle", "isolated; crossed by the -20% line"),
    (765, 468): ("open triangle", "isolated"),
}
# Touching pairs split by hand (centre read on the 3x zoom, native px).
SPLIT = {
    (768, 390): [("open inverted triangle", 772.7, 383.3), ("open inverted triangle", 765.0, 394.3)],
    (554, 654): [("open circle", 550.7, 656.0), ("open circle", 557.3, 650.0)],
}
TRIANGLES = {"open triangle": +1, "filled triangle": +1,
             "open inverted triangle": -1, "filled inverted triangle": -1}


def page_image(pdf: str) -> np.ndarray:
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["pdfimages", "-png", "-f", str(PDF_PAGE), "-l", str(PDF_PAGE), pdf,
                        f"{tmp}/p"], check=True)
        png = sorted(Path(tmp).glob("p-*.png"))[0]
        return np.array(Image.open(png).convert("L")) < 128


def calibrate(a):
    """Frame, ticks and the two affine maps log10(value) <- px; returns (fx, fy, residuals)."""
    col = 240 + int(np.argmax(a[220:950, 240:251].sum(0)))
    w = a[220:950, col - 2:col + 3].sum(0)
    x_frame = col - 2 + float((np.arange(5) * w).sum() / w.sum())
    row = 955 + int(np.argmax(a[955:968, 260:1000].sum(1)))
    w = a[row - 2:row + 3, 260:1000].sum(1)
    y_frame = row - 2 + float((np.arange(5) * w).sum() / w.sum())
    xt = np.where(a[940:958, 250:1010].sum(0) >= 8)[0] + 250
    yt = np.where(a[200:955, 247:265].sum(1) >= 8)[0] + 200
    xt = [g.mean() for g in np.split(xt, np.where(np.diff(xt) > 1)[0] + 1)]
    yt = [g.mean() for g in np.split(yt, np.where(np.diff(yt) > 1)[0] + 1)]
    if len(xt) != 2 or len(yt) != 2:
        raise RuntimeError(f"ticks not found: {xt}, {yt}")
    tx = np.array([[x_frame, -2], [xt[0], -1], [xt[1], 0]])
    ty = np.array([[y_frame, -2], [yt[1], -1], [yt[0], 0]])
    fx = np.polyfit(tx[:, 0], tx[:, 1], 1)
    fy = np.polyfit(ty[:, 0], ty[:, 1], 1)
    res_x = (tx[:, 1] - np.polyval(fx, tx[:, 0])) / fx[0]
    res_y = (ty[:, 1] - np.polyval(fy, ty[:, 0])) / fy[0]
    return fx, fy, tx, ty, res_x, res_y


def main(pdf):
    a = page_image(pdf)
    fx, fy, tx, ty, res_x, res_y = calibrate(a)
    r0, r1 = ROWS
    c0, c1 = COLS
    sub = a[r0:r1, c0:c1]
    R, C = np.mgrid[r0:r1, c0:c1]
    D = np.polyval(fy, R) - np.polyval(fx, C)
    lines = {}
    for k, (lo, hi) in {"solid": (-0.02, 0.03), "plus": (0.07, 0.11), "minus": (-0.115, -0.075)}.items():
        cen = D[sub & (D > lo) & (D < hi)].mean()
        for w in (0.012, 0.008, 0.006):
            cen = D[sub & (np.abs(D - cen) < w)].mean()
        lines[k] = float(cen)
    keep = sub.copy()
    for v in lines.values():
        keep &= ~(np.abs(D - v) < LINE_HALF)
    lab, _ = ndi.label(ndi.binary_dilation(keep, iterations=2))
    rows = []
    for i, sl in enumerate(ndi.find_objects(lab)):
        m = (lab[sl] == i + 1) & keep[sl]
        if m.sum() < MIN_AREA:
            continue
        rr, cc = np.nonzero(m)
        rr = rr + sl[0].start + r0
        cc = cc + sl[1].start + c0
        x = (cc.min() + cc.max()) / 2
        y = (rr.min() + rr.max()) / 2
        h = rr.max() - rr.min() + 1
        key = min(list(HAND) + list(SPLIT), key=lambda p: (p[0] - x)**2 + (p[1] - y)**2)
        near = (key[0] - x)**2 + (key[1] - y)**2 < 6**2
        dpx = np.polyval(fy, rr) - np.polyval(fx, cc)
        if near and key in HAND:
            sym, flag = HAND[key]
            y += TRIANGLES.get(sym, 0) * h / 6.0
            rows.append(("marker", sym, flag, x, y, np.nan, np.nan))
        elif near and key in SPLIT:
            for sym, xs, ys in SPLIT[key]:
                rows.append(("marker", sym, "touching pair; split by hand", xs, ys, np.nan, np.nan))
        elif m.sum() > 150 or 500 < x < 800:
            rows.append(("cluster", "several", f"{m.sum()} px of touching markers" if m.sum() > 150 else f"{m.sum()} px fragment of a marker chain", x, y,
                         float(dpx.min()), float(dpx.max())))
    with open(OUT, "w") as fh:
        fh.write("# Saez and Carbonell (1985) Figure 8 (p. 60): liquid holdup, calculated against "
                 "experimental, log-log. Digitised by data/digitise.py; see the sidecar.\n")
        fh.write("kind,symbol,flag,x_px,y_px,exp,cal,offset_min_dec,offset_max_dec\n")
        for kind, sym, flag, x, y, dmin, dmax in sorted(rows, key=lambda t: t[3]):
            ex, ca = 10**np.polyval(fx, x), 10**np.polyval(fy, y)
            if kind == "cluster":
                fh.write(f"{kind},{sym},{flag},{x:.1f},{y:.1f},,,{dmin:.5f},{dmax:.5f}\n")
            else:
                fh.write(f"{kind},{sym},{flag},{x:.1f},{y:.1f},{ex:.5f},{ca:.5f},,\n")
        for k, v in lines.items():
            fh.write(f"line,{k},drawn line offset from the tick calibration,,,,,{v:.5f},{v:.5f}\n")
    print("x ticks (px, log10):", tx.tolist(), "residual px", np.round(res_x, 2).tolist())
    print("y ticks (px, log10):", ty.tolist(), "residual px", np.round(res_y, 2).tolist())
    print("px per decade x %.2f, y %.2f" % (1 / fx[0], -1 / fy[0]))
    print("lines (decades):", {k: round(v, 5) for k, v in lines.items()})
    print(OUT.read_text())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
