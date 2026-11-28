"""Digitised readings of Al-Dahhan and Dudukovic (1995), Figures 4 to 7, and their calibration.

The readings are marker centres (centres of the marker's bounding box) and line
points in page pixels of the PDF pages rendered at 300 ppi (`pdftoppm -r 300`), the
native resolution of the scan. Markers and their boxes were found as connected
components of the thresholded image; touching markers were separated by eye on
pixel printouts, rings of Fig. 6 by a ring-template score. Each axis is calibrated
by a least-squares fit of the value (linear axes) or of its logarithm (Figs 4a, 5a
and 7, whose tick spacing is logarithmic although the labels read linear) on the
pixel coordinate of the printed ticks. Run from the page directory:

    python data/digitise.py

It writes al-dahhan-1995-fig4.csv to -fig7.csv and prints the tick-fit residuals
that the sidecars quote. No page image is stored in the repository.
"""
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SRC = "Al-Dahhan and Dudukovic (1995)"


class Axis:
    """Least-squares map from a pixel coordinate to a value, linear or logarithmic."""

    def __init__(self, ticks, log):
        v, p = np.array(ticks, dtype=float).T
        self.log = log
        y = np.log(v) if log else v
        self.b, self.a = np.polyfit(p, y, 1)
        self.res_px = (y - (self.a + self.b * p)) / self.b     # residual in pixels
        self.px_per_unit = 1.0 / abs(self.b)                   # px per unit (ln unit if log)

    def __call__(self, px):
        y = self.a + self.b * np.asarray(px, dtype=float)
        return np.exp(y) if self.log else y


# ---- Figure 4 (journal p. 2385, PDF page 9): bed of spheres.
# (a) eta_CE against L, log-log; ticks: frame corners 0.5 and 3.0 (x), 1.0 and 0.5 (y).
F4A_X = Axis([(0.5, 415.5), (1.0, 658), (2.0, 901), (3.0, 1044.5)], log=True)
F4A_Y = Axis([(1.0, 2199), (0.9, 2295), (0.8, 2400.5), (0.7, 2521.5), (0.6, 2660.5), (0.5, 2823)], log=True)
# (b) (dP/Z)/(rho_L g) against L, linear.
F4B_X = Axis([(0, 1327), (0.5, 1431), (1.0, 1533.5), (1.5, 1635), (2.0, 1737.5), (2.5, 1840.5), (3.0, 1942.5)], log=False)
F4B_Y = Axis([(0, 1565), (2, 1411), (4, 1257), (6, 1103), (8, 949)], log=False)
# (case, shape, px_x, px_y, reading SD in px, flag)
F4A = [
    (2, "star", 500.5, 2631.5, 1, ""), (2, "star", 704, 2427.5, 1, ""), (2, "star", 821, 2311, 1, ""),
    (2, "star", 857, 2279, 1, ""), (2, "star", 1002, 2198.5, 2, "frame"),
    (3, "cross", 615, 2487.5, 1, ""), (3, "cross", 700, 2388, 1, ""), (3, "cross", 796.5, 2330, 2.5, "overlap"),
    (3, "cross", 876.5, 2233, 1, ""),
    (4, "open square", 557, 2566.5, 1, ""), (4, "open square", 729.5, 2396.5, 1, ""),
    (4, "open square", 802.5, 2333, 2.5, "overlap"), (4, "open square", 909.5, 2247.5, 1, ""),
    (5, "filled square", 583, 2444.5, 1, ""), (5, "filled square", 742.5, 2312, 1, ""),
    (5, "filled square", 794, 2265, 1, ""), (5, "filled square", 888.5, 2198.5, 2, "frame"),
]
F4B = [
    (1, "filled triangle", 1412, 1565, 2.5, "axis"), (1, "filled triangle", 1502.5, 1565, 2.5, "axis"),
    (1, "filled triangle", 1555.5, 1565, 2.5, "axis"), (1, "filled triangle", 1644.5, 1565, 2.5, "axis"),
    (1, "filled triangle", 1710.5, 1565, 2.5, "axis"),
    (2, "star", 1459.5, 1538, 1, ""), (2, "star", 1562, 1536.5, 2.5, "overlap"), (2, "star", 1654, 1531, 2.5, "overlap"),
    (2, "star", 1689, 1523, 1, ""), (2, "star", 1873.5, 1516, 1, ""),
    (3, "cross", 1509, 1487.5, 1, ""), (3, "cross", 1558.5, 1473, 1, ""), (3, "cross", 1631.5, 1468.5, 1, ""),
    (3, "cross", 1709.5, 1462, 1, ""),
    (4, "open square", 1482, 1536.5, 1, ""), (4, "open square", 1579, 1530.5, 2.5, "overlap"),
    (4, "open square", 1637.5, 1529, 2.5, "overlap"), (4, "open square", 1747, 1527, 1, ""),
    (5, "filled square", 1493, 1291.5, 1, ""), (5, "filled square", 1587.5, 1249, 1, ""),
    (5, "filled square", 1629.5, 1208.5, 1, ""), (5, "filled square", 1723, 1195, 1, ""),
]

# ---- Figure 5 (journal p. 2386, PDF page 10): bed of extrudates. Case "5b" is the
# 4.93 MPa series (open triangle), printed as a second "5" in the legend.
F5A_X = Axis([(0.5, 445), (0.6, 515.5), (1.0, 713), (2.0, 981.5)], log=True)
F5A_Y = Axis([(1.0, 230.5), (0.9, 325), (0.8, 431.5), (0.7, 553.5), (0.6, 693.5), (0.5, 857)], log=True)
F5B_X = Axis([(0, 454.5), (0.5, 555.5), (1.0, 657.5), (1.5, 758.5), (2.0, 860.5), (2.5, 962.5), (3.0, 1066.5)], log=False)
F5B_Y = Axis([(0, 1752), (2, 1576), (4, 1400.5), (6, 1226)], log=False)
F5A = [
    (2, "star", 596, 658.5, 1, ""), (2, "star", 762, 498.5, 1, ""), (2, "star", 854, 443, 1, ""), (2, "star", 935, 365.5, 1, ""),
    (3, "cross", 641, 575, 1, ""), (3, "cross", 744, 485, 1, ""), (3, "cross", 821, 414.5, 1, ""),
    (3, "cross", 939.5, 320.5, 2.5, "overlap"),
    (4, "open square", 650, 594, 1, ""), (4, "open square", 792.5, 468, 1, ""), (4, "open square", 847.5, 423.5, 1, ""),
    (4, "open square", 927, 335, 2.5, "overlap"),
    (5, "filled square", 669, 466.5, 1, ""), (5, "filled square", 757.5, 384, 1, ""),
    (5, "filled square", 927, 251.5, 1, ""), (5, "filled square", 990, 230.5, 2, "frame"),
    ("5b", "open triangle", 551, 610, 1, ""), ("5b", "open triangle", 825, 346.5, 1, ""),
    ("5b", "open triangle", 948.5, 242, 1, ""), ("5b", "open triangle", 1038.5, 230, 2, "frame"),
]
F5B = [
    (2, "star", 604, 1746, 2.5, "axis"), (2, "star", 685.5, 1745.5, 2.5, "axis"),
    (2, "star", 747, 1746, 3, "overlap; axis"), (2, "star", 815, 1744, 3, "overlap; axis"),
    (3, "cross", 623, 1692.5, 1, ""), (3, "cross", 675, 1677, 1, ""), (3, "cross", 724, 1666, 1, ""),
    (3, "cross", 819.5, 1654, 1, ""),
    (4, "open square", 626.5, 1736, 1, ""), (4, "open square", 704, 1731, 1, ""),
    (4, "open square", 741.5, 1728.5, 2.5, "overlap"), (4, "open square", 807.5, 1726, 2.5, "overlap"),
    (5, "filled square", 636, 1437.5, 1, ""), (5, "filled square", 683, 1409.5, 1, ""),
    (5, "filled square", 808, 1377, 1, ""), (5, "filled square", 870, 1367, 1, ""),
    ("5b", "open triangle", 589, 1475, 1, ""), ("5b", "open triangle", 727.5, 1419, 1, ""),
    ("5b", "open triangle", 828.5, 1418.5, 1, ""), ("5b", "open triangle", 925.5, 1406, 1, ""),
]

# ---- Figure 6 (journal p. 2387, PDF page 11): eta_CE experimental against predicted,
# linear, open rings. Calibrated on the labelled inner ticks 0.2 to 0.8.
F6_X = Axis([(0.2, 526), (0.4, 654.5), (0.6, 783), (0.8, 913)], log=False)
F6_Y = Axis([(0.8, 370), (0.6, 500), (0.4, 629.5), (0.2, 758.5)], log=False)
F6 = [
    (777.0, 492.6, 1, ""), (800.5, 455.8, 1, ""), (806.6, 483.5, 1, ""), (831.0, 454.7, 1, ""),
    (840.6, 416.5, 1, ""), (844.4, 445.0, 1, ""), (855.9, 461.5, 1, ""), (859.1, 383.5, 1.5, "overlap"),
    (869.5, 392.9, 1.5, "overlap"), (870.0, 406.8, 1.5, "overlap"), (873.8, 367.3, 1.5, "overlap"),
    (884.8, 360.1, 1.5, "overlap"), (890.0, 388.5, 2, "overlap; double"), (891.5, 400.5, 2, "overlap; double"),
    (907.5, 374.5, 2, "overlap"), (909.5, 326.4, 1.5, "overlap"), (909.5, 388.3, 1.5, "overlap"),
    (916.2, 363.7, 1.5, "overlap"), (918.0, 314.1, 1.5, "overlap"), (930.6, 359.1, 1, ""),
    (932.1, 326.0, 1.5, "overlap"), (947.5, 329.2, 1.5, "overlap"),
    (953.6, 312.9, 3, "double ring; two coincident markers"), (953.6, 312.9, 3, "double ring; two coincident markers"),
    (955.7, 340.8, 1.5, "overlap"), (966.5, 270.8, 1.5, "overlap"), (976.8, 262.3, 1.5, "overlap"),
    (987.9, 317.5, 1, ""), (992.0, 301.5, 1, ""),
    (1020.5, 238.5, 2, "frame"),
    (1040.5, 238.0, 3, "corner stack on both frames; count of 3 or 4 not decidable"),
    (1040.5, 245.5, 3, "corner stack on both frames; count of 3 or 4 not decidable"),
    (1040.5, 254.5, 3, "corner stack on both frames; count of 3 or 4 not decidable"),
]

# ---- Figure 7 (journal p. 2387, PDF page 11): extrudates, eta_CE against L, log-log.
F7_X = Axis([(0.5, 422.5), (1.0, 670.5), (2.0, 919.5), (3.0, 1065)], log=True)
F7_Y = Axis([(1.0, 1274.5), (0.9, 1407), (0.8, 1554), (0.7, 1722), (0.6, 1913.5)], log=True)
F7 = [  # (kind, series, shape, px_x, px_y, SD px, flag)
    ("marker", "exp 0.31 MPa", "open circle", 564.5, 1867.0, 1, ""),
    ("marker", "exp 0.31 MPa", "open circle", 716.0, 1647.5, 1.5, "on the solid line"),
    ("marker", "exp 0.31 MPa", "open circle", 800.0, 1571.5, 1.5, "on the solid line"),
    ("marker", "exp 0.31 MPa", "open circle", 875.3, 1463.8, 1.5, "on the solid line"),
    ("marker", "exp 3.55 MPa", "open square", 629.5, 1600.5, 1, ""),
    ("marker", "exp 3.55 MPa", "open square", 711.5, 1490.0, 1, ""),
    ("marker", "exp 3.55 MPa", "open square", 866.5, 1305.5, 1, ""),
] + [("line", "pred 0.31 MPa", "solid", x, y, 0.7, "line centre row") for x, y in (
    (577, 1817.5), (597, 1794.5), (617, 1771.0), (637, 1748.0), (657, 1724.0), (677, 1700.5), (697, 1676.0),
    (740, 1624.0), (760, 1600.5), (780, 1577.0), (823, 1526.0), (843, 1503.0), (863, 1479.5))] + \
    [("line", "pred 3.55 MPa", "dashed", x, y, 1.0, "dash centroid") for x, y in (
        (643.1, 1541.2), (687.6, 1487.2), (724.1, 1443.8), (765.1, 1397.0), (801.3, 1351.7), (849.1, 1298.2))]


def write(name, header, cols, rows):
    with open(HERE / name, "w") as f:
        f.write(f"# {header}\n" + ",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(v) for v in r) + "\n")


def sd_val(ax, px, sd_px):
    """Reading SD in the value (relative for a log axis, absolute for a linear one)."""
    return sd_px / ax.px_per_unit


def panel_rows(panel, pts, ax_x, ax_y, yname):
    out = []
    for case, shape, px, py, sd, flag in pts:
        L, y = float(ax_x(px)), float(ax_y(py))
        out.append((panel, case, shape, f"{L:.5f}", yname, f"{y:.5f}", px, py,
                    f"{sd_val(ax_x, px, sd):.5f}", f"{sd_val(ax_y, py, sd):.5f}", flag))
    return out


COLS = ["panel", "case", "shape", "L_kg_m2s", "quantity", "value", "px_x", "px_y", "sd_L", "sd_value", "flag"]


def main():
    for name, ax in [("4a x", F4A_X), ("4a y", F4A_Y), ("4b x", F4B_X), ("4b y", F4B_Y), ("5a x", F5A_X),
                     ("5a y", F5A_Y), ("5b x", F5B_X), ("5b y", F5B_Y), ("6 x", F6_X), ("6 y", F6_Y),
                     ("7 x", F7_X), ("7 y", F7_Y)]:
        print(f"Fig. {name}: {ax.px_per_unit:8.2f} px per {'ln unit' if ax.log else 'unit'}, "
              f"largest tick residual {np.abs(ax.res_px).max():.2f} px")
    write("al-dahhan-1995-fig4.csv",
          f"{SRC} Fig. 4 (p. 2385), bed of spheres: (a) eta_CE and (b) (dP/Z)/(rho_L g) against L. "
          "Digitised by data/digitise.py; see the sidecar. sd_L and sd_value: relative on log axes (a), absolute on linear (b).",
          COLS, panel_rows("a", F4A, F4A_X, F4A_Y, "eta_CE") + panel_rows("b", F4B, F4B_X, F4B_Y, "dpz_rhog"))
    write("al-dahhan-1995-fig5.csv",
          f"{SRC} Fig. 5 (p. 2386), bed of extrudates: (a) eta_CE and (b) (dP/Z)/(rho_L g) against L. "
          "Case 5b is the 4.93 MPa series (open triangle). Digitised by data/digitise.py; see the sidecar.",
          COLS, panel_rows("a", F5A, F5A_X, F5A_Y, "eta_CE") + panel_rows("b", F5B, F5B_X, F5B_Y, "dpz_rhog"))
    rows6 = [(i, f"{float(F6_X(px)):.5f}", f"{float(F6_Y(py)):.5f}", px, py,
              f"{sd / F6_X.px_per_unit:.5f}", f"{sd / F6_Y.px_per_unit:.5f}", flag)
             for i, (px, py, sd, flag) in enumerate(F6)]
    write("al-dahhan-1995-fig6.csv",
          f"{SRC} Fig. 6 (p. 2387): eta_CE experimental against eta_CE predicted, both beds, open rings. "
          "Digitised by data/digitise.py; see the sidecar.",
          ["ring", "eta_pred", "eta_exp", "px_x", "px_y", "sd_pred", "sd_exp", "flag"], rows6)
    rows7 = [(k, s, sh, f"{float(F7_X(px)):.5f}", f"{float(F7_Y(py)):.5f}", px, py,
              f"{sd / F7_X.px_per_unit:.5f}", f"{sd / F7_Y.px_per_unit:.5f}", flag)
             for k, s, sh, px, py, sd, flag in F7]
    write("al-dahhan-1995-fig7.csv",
          f"{SRC} Fig. 7 (p. 2387): extrudates, eta_CE against L, measured markers and the two predicted lines. "
          "Log-log axes. Digitised by data/digitise.py; see the sidecar. sd_L and sd_eta relative.",
          ["kind", "series", "shape", "L_kg_m2s", "eta_CE", "px_x", "px_y", "sd_L", "sd_eta", "flag"], rows7)


if __name__ == "__main__":
    main()
