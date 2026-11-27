"""Digitised readings of Zahradnik et al. (1997), Figures 3(a), 11(a) and 12(b), and their calibration.

The readings below are marker centres and drawn-line centres in pixels of the native
page bitmaps of the publisher PDF (`pdfimages -png`, 2280 x 3255 px, 300 ppi, 1-bit),
PDF pages 5 (Fig. 3a), 11 (Fig. 11a) and 12 (Fig. 12b). Markers were found by a
ring template (open circles, radius 5.5 to 9 px) or a disk template (filled circles,
radius 7.5 px) and checked by eye against the page; every marker on the page was
found and no legend glyph was kept. A centre is the centroid of the dark pixels within
the marker radius, so a drawn line through a marker pulls it by about a pixel.

Each scan is rotated by about half a degree. The frame lines give the rotation: an
iso-abscissa line runs parallel to the left frame, an iso-ordinate line parallel to the
bottom frame (slopes below, fitted to the frame centroids). Each axis is then calibrated
by a straight-line fit of the printed ticks. Run from the page directory:

    python data/digitise.py

It writes the three figure CSVs and prints the tick-fit residuals the sidecars quote.
No page image is stored in the repository.
"""
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent

# ---- Figure 11(a) (journal p. 3821, PDF page 11): u_s/u_binf against eps_G.
# Left frame: col = -0.01018 row + 371.32; bottom frame: row = 0.00979 col + 2760.31.
F11_ROT = (0.01018, 2755.0, 0.00979, 355.0)   # (left slope, tick row, bottom slope, tick col)
F11_XT = [(0.1, 539.5), (0.2, 733.5), (0.3, 927.5)]          # bottom-frame ticks (u after rotation)
F11_YT = [(0.0, 2763.79), (1.0, 2515.54), (2.0, 2266.50)]    # bottom frame, right tick 1, left tick 2 (v)
F11_MARKERS = [(384.2, 2522.6), (433.1, 2530.3), (483.5, 2535.7), (582.1, 2548.2), (647.0, 2555.7),
               (701.1, 2563.2), (754.3, 2569.2), (790.4, 2562.7), (841.3, 2560.4), (922.7, 2555.9),
               (983.6, 2545.2), (1020.2, 2527.2), (1029.8, 2515.9), (1026.0, 2498.4), (1005.2, 2473.8),
               (966.0, 2439.1), (892.4, 2353.0), (850.1, 2275.9), (838.6, 2225.9), (860.7, 2205.2),
               (884.3, 2191.2)]
# Drawn Table 2 lines 1 to 5: centre of the line's dark run at eps 0.05 to 0.30, where no
# marker covers it and no other line lies within 0.03 in u_s/u_binf (line number, px, py).
F11_LINES = [(1, 444.9, 2517.1), (2, 445.0, 2504.5), (1, 493.4, 2516.1), (2, 493.6, 2497.4),
             (1, 541.9, 2518.1), (2, 542.2, 2490.4), (4, 541.6, 2550.4), (1, 590.4, 2517.0),
             (2, 590.8, 2482.0), (1, 638.9, 2518.5), (2, 639.4, 2472.6), (4, 638.4, 2566.8),
             (1, 687.4, 2520.0), (2, 688.0, 2463.2), (4, 686.8, 2574.8), (5, 687.1, 2548.6),
             (1, 735.9, 2519.5), (2, 736.6, 2453.3), (4, 735.3, 2582.3), (5, 735.6, 2552.6),
             (1, 784.4, 2520.9), (2, 785.2, 2444.2), (3, 783.8, 2574.4), (4, 783.7, 2590.6),
             (1, 832.9, 2520.4), (2, 833.8, 2433.0), (3, 832.3, 2580.6), (4, 832.1, 2598.9),
             (1, 881.4, 2521.4), (2, 882.4, 2421.6), (3, 880.7, 2589.3), (4, 880.5, 2609.0),
             (5, 881.0, 2562.0), (1, 929.9, 2521.4), (2, 931.0, 2409.8), (3, 929.1, 2599.4),
             (4, 928.9, 2620.0), (5, 929.4, 2566.4)]

# ---- Figure 3(a) (journal p. 3815, PDF page 5): eps_G against u_OG, plate A open circles,
# plate B filled circles (the k_L a triangles are not read).
# Left frame: col = 0.01254 row + 311.11; bottom frame: row = -0.00815 col + 1093.14.
F3_ROT = (-0.01254, 1082.0, -0.00815, 335.0)
F3_XT = [(0.04, 506.5), (0.05, 551.5), (0.10, 780.5), (0.125, 895.5), (0.15, 1010.5)]
F3_YT = [(0.0, 1090.4), (0.1, 861.5), (0.2, 627.5), (0.3, 397.5)]
F3_OPEN = [(396.9, 910.3), (447.6, 760.9), (485.8, 643.1), (522.3, 550.3), (561.5, 440.0),
           (596.8, 378.2), (642.7, 323.1), (689.7, 271.4), (733.4, 281.0), (762.9, 329.2),
           (811.5, 422.8), (847.3, 478.8), (887.3, 497.6), (979.9, 452.7)]
F3_FILLED = [(404.7, 945.3), (466.4, 862.1), (565.1, 748.3), (646.2, 677.3), (732.7, 615.7),
             (808.5, 555.2), (889.0, 509.0), (976.1, 459.2)]

# ---- Figure 12(b) (journal p. 3822, PDF page 12): u_OG/eps_G against u_OG, plate B.
# Left frame: col = 0.00480 row + 1215.25; bottom frame: row = -0.000526 col + 2760.91.
F12_ROT = (-0.0048, 2750.0, -0.000526, 1235.0)
F12_XT = [(0.05, 1456.5), (0.10, 1686.5), (0.15, 1915.5)]
F12_YT = [(0.0, 2760.3), (0.1, 2666.0), (0.2, 2571.5), (0.3, 2480.5), (0.4, 2387.0), (0.5, 2294.5),
          (0.6, 2200.5), (0.7, 2106.0), (0.8, 2012.5)]
F12_MARKERS = [(1310.0, 2489.2), (1367.7, 2465.7), (1474.0, 2417.8), (1556.5, 2380.2), (1639.2, 2345.0),
               (1720.9, 2315.5), (1803.7, 2282.6), (1886.2, 2250.1)]


def calibration(rot, xt, yt):
    """Map (px, py) to axis values: rotate by the frame slopes, then the tick fits."""
    sl, row0, sb, col0 = rot
    xv, xp = np.array(xt).T
    yv, yp = np.array(yt).T
    cx, cy = np.polyfit(xv, xp, 1), np.polyfit(yv, yp, 1)
    rx, ry = xp - np.polyval(cx, xv), yp - np.polyval(cy, yv)

    def to_axes(px, py):
        u = np.asarray(px) + sl * (np.asarray(py) - row0)
        v = np.asarray(py) - sb * (np.asarray(px) - col0)
        return (u - cx[1]) / cx[0], (v - cy[1]) / cy[0]
    return to_axes, cx[0], cy[0], np.abs(rx).max(), np.abs(ry).max()


def write(name, header, cols, rows):
    with open(HERE / name, "w") as f:
        f.write(header + "\n" + ",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(r) + "\n")


if __name__ == "__main__":
    t11, sx, sy, rx, ry = calibration(F11_ROT, F11_XT, F11_YT)
    print(f"Fig. 11a: {sx:.1f} px per unit eps, {-sy:.2f} px per unit y; tick residuals {rx:.2f}, {ry:.2f} px")
    rows = []
    for i, (px, py) in enumerate(F11_MARKERS, 1):
        e, y = t11(px, py)
        rows.append(["marker", str(i), f"{e:.4f}", f"{y:.4f}", f"{px:.1f}", f"{py:.1f}"])
    for k, px, py in F11_LINES:
        e, y = t11(px, py)
        rows.append(["line", str(k), f"{e:.4f}", f"{y:.4f}", f"{px:.1f}", f"{py:.1f}"])
    write("zahradnik-1997-fig11a.csv",
          "# Zahradnik et al. (1997) Fig. 11(a) (p. 3821): slip velocity ratio against gas holdup, air-water, "
          "D_r 0.15 m, H0/D_r 10, plate A. Digitised by data/digitise.py; see the sidecar.",
          ["kind", "series", "eps_G", "us_over_ubinf", "px", "py"], rows)

    t3, sx, sy, rx, ry = calibration(F3_ROT, F3_XT, F3_YT)
    print(f"Fig. 3a: {sx:.1f} px per m/s, {-sy:.1f} px per unit eps; tick residuals {rx:.2f}, {ry:.2f} px")
    rows = []
    for plate, pts in (("A", F3_OPEN), ("B", F3_FILLED)):
        for px, py in pts:
            u, e = t3(px, py)
            rows.append([plate, f"{u:.5f}", f"{e:.4f}", f"{px:.1f}", f"{py:.1f}"])
    write("zahradnik-1997-fig3a.csv",
          "# Zahradnik et al. (1997) Fig. 3(a) (p. 3815): gas holdup against superficial gas velocity, air-water, "
          "D_r 0.15 m, H0/D_r 10, plates A (open) and B (filled). Digitised by data/digitise.py; see the sidecar.",
          ["plate", "u_G", "eps_G", "px", "py"], rows)

    t12, sx, sy, rx, ry = calibration(F12_ROT, F12_XT, F12_YT)
    print(f"Fig. 12b: {sx:.1f} px per m/s, {-sy:.1f} px per m/s; tick residuals {rx:.2f}, {ry:.2f} px")
    rows = []
    for px, py in F12_MARKERS:
        u, y = t12(px, py)
        rows.append([f"{u:.5f}", f"{y:.4f}", f"{px:.1f}", f"{py:.1f}"])
    write("zahradnik-1997-fig12b.csv",
          "# Zahradnik et al. (1997) Fig. 12(b) (p. 3822): drift-flux plot, bubble column D_r 0.15 m, H0/D_r 10, "
          "plate B, semibatch. Digitised by data/digitise.py; see the sidecar.",
          ["u_G", "uG_over_eps", "px", "py"], rows)
