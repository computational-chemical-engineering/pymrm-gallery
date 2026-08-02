"""Figure 3 (journal page 1383): (Sh-2)/Sc^(1/3) vs Re for the liquid-phase data.

This is the figure alpha and beta were FITTED on.  Its ordinate is already the
reduced group, so a free refit of Sh = 2 + alpha Sc^(1/3) Re^beta against it is a
plain power law y = alpha Re^beta - no Schmidt number is needed and none is
printed.
"""
import os

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

S = os.environ.get("A34_RENDERS", "renders")
PATH = f"{S}/fig3-full.png"

# decade positions measured on the printed ticks (see calibrate() below).
# The ordinate calibration is COLUMN-DEPENDENT: the render is very slightly
# rotated, so the left-hand decade ticks sit ~5 px above the right-hand ones and
# a single row origin is wrong by up to half that at either end of the plot.
COL0 = PXDEC_X = None
ROW0 = PXDEC_Y = None                       # the flat (both-sides averaged) fit
ROW0_L = ROW0_R = PXDEC_Y_L = PXDEC_Y_R = None
TICKCOL_L = TICKCOL_R = None                # columns the two tick strips sample
FRAME = dict(left=142.5, right=1753.5, top=135.5, bottom=1148.5)

TEXT_BOXES = [
    (170, 505, 220, 840),      # legend, left column
    (170, 378, 975, 1665),     # legend, right column
    (915, 1020, 875, 1640),    # "(Sh - 2)/Sc^(1/3) = 1 1 Re^0 6"
]


def load():
    return 1.0 - np.asarray(Image.open(PATH).convert("L"), float) / 255.0


def _peaks(profile, min_sep, thr):
    idx = np.where(profile > thr)[0]
    if len(idx) == 0:
        return []
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        (cur if i - cur[-1] <= min_sep else groups.append(cur) or cur.clear() or cur).append(i) \
            if False else None
        if i - cur[-1] <= min_sep:
            cur.append(i)
        else:
            groups.append(cur); cur = [i]
    groups.append(cur)
    return [float(np.sum(np.array(g) * profile[g]) / profile[g].sum()) for g in groups]


def calibrate(ink, verbose=False):
    """Least-squares decade scale from the printed major ticks plus the frame.

    x: the frame's left and right edges are the 1 and 10^4 gridlines; the three
       interior decades are read off the ticks hanging from the top axis, which
       is the one edge no data touches.
    y: the left and right axes each carry all four decade ticks (y = 1, 10, 100,
       1000), and the two sides DISAGREE - the left ticks sit about 5 px above
       the right ones, because the page render is rotated by ~0.2 degrees.  Each
       side is therefore fitted separately and the ordinate origin and decade
       height are linearly interpolated between the two tick columns.  Averaging
       the two sides instead (what the 2026-08-02 pass did) is right in the
       middle of the plot and wrong by ~2.6 px at either end, which is worth
       -0.005 in beta and +2.5 % in alpha; see the drawn-line control on the
       page.
    """
    global COL0, PXDEC_X, ROW0, PXDEC_Y
    global ROW0_L, ROW0_R, PXDEC_Y_L, PXDEC_Y_R, TICKCOL_L, TICKCOL_R
    top = ink[int(FRAME["top"]) + 4:int(FRAME["top"]) + 19, :] > 0.5
    xt = _peaks(top.sum(0).astype(float), 6, 8)
    xd, xc = [0.0, 4.0], [FRAME["left"], FRAME["right"]]
    for k, want in ((1, FRAME["left"] + 402.7), (2, FRAME["left"] + 805.4),
                    (3, FRAME["left"] + 1208.1)):
        cand = [v for v in xt if abs(v - want) < 12]
        if cand:
            xd.append(float(k)); xc.append(min(cand, key=lambda v: abs(v - want)))
    A = np.polyfit(xd, xc, 1)
    PXDEC_X, COL0 = float(A[0]), float(A[1])

    lo_l, hi_l = int(FRAME["left"]) + 4, int(FRAME["left"]) + 19
    lo_r, hi_r = int(FRAME["right"]) - 18, int(FRAME["right"]) - 3
    TICKCOL_L, TICKCOL_R = 0.5 * (lo_l + hi_l - 1), 0.5 * (lo_r + hi_r - 1)
    yl = _peaks((ink[:, lo_l:hi_l] > 0.5).sum(1).astype(float), 6, 8)
    yr = _peaks((ink[:, lo_r:hi_r] > 0.5).sum(1).astype(float), 6, 8)
    # the four decade ticks, wanted at the nominal positions of y = 1 .. 1000
    wants = [FRAME["bottom"] - 4.0, FRAME["bottom"] - 341.7,
             FRAME["bottom"] - 679.4, FRAME["top"] + 0.5]
    side = []
    for peaks in (yl, yr):
        kk, rr = [], []
        for k, want in enumerate(wants):
            cand = [v for v in peaks if abs(v - want) < 12]
            if cand:
                kk.append(float(k)); rr.append(min(cand, key=lambda v: abs(v - want)))
        B = np.polyfit(kk, rr, 1)
        side.append((float(B[1]), float(-B[0]), kk, rr,
                     float(np.std(np.asarray(rr) - np.polyval(B, kk)))))
    (ROW0_L, PXDEC_Y_L, kl, rl, resl), (ROW0_R, PXDEC_Y_R, kr, rr_, resr) = side
    ROW0 = 0.5 * (ROW0_L + ROW0_R)
    PXDEC_Y = 0.5 * (PXDEC_Y_L + PXDEC_Y_R)
    if verbose:
        print(f"  x: decades {xd} at columns {np.round(xc,1)}")
        print(f"     Re = 1 at column {COL0:.2f}, {PXDEC_X:.2f} px per decade")
        print(f"  y left  (col {TICKCOL_L:.0f}): decades {kl} at rows "
              f"{np.round(rl,2)} -> y=1 at {ROW0_L:.2f}, {PXDEC_Y_L:.2f} px/dec "
              f"(resid {resl:.2f} px)")
        print(f"  y right (col {TICKCOL_R:.0f}): decades {kr} at rows "
              f"{np.round(rr_,2)} -> y=1 at {ROW0_R:.2f}, {PXDEC_Y_R:.2f} px/dec "
              f"(resid {resr:.2f} px)")
        print(f"     skew: the y = 1 line drops {ROW0_R - ROW0_L:+.2f} px from the "
              f"left tick column to the right one")
    return COL0, PXDEC_X, ROW0, PXDEC_Y


def re_of_col(c):
    return 10.0 ** ((np.asarray(c, float) - COL0) / PXDEC_X)


def col_of_re(re):
    return COL0 + PXDEC_X * np.log10(np.asarray(re, float))


def _yscale(c):
    """Ordinate origin and decade height at column `c`, interpolated between the
    two tick columns.  Both arguments of the ordinate calibration are therefore
    functions of the column: that is the skew correction."""
    t = (np.asarray(c, float) - TICKCOL_L) / (TICKCOL_R - TICKCOL_L)
    return (ROW0_L + (ROW0_R - ROW0_L) * t,
            PXDEC_Y_L + (PXDEC_Y_R - PXDEC_Y_L) * t)


def y_of_row(r, c):
    r0, py = _yscale(c)
    return 10.0 ** ((r0 - np.asarray(r, float)) / py)


def row_of_y(y, c):
    r0, py = _yscale(c)
    return r0 - py * np.log10(np.asarray(y, float))


def y_of_row_flat(r):
    """The 2026-08-02 calibration: one row origin for the whole plot.  Kept so
    the page can show what the skew correction moved."""
    return 10.0 ** ((ROW0 - np.asarray(r, float)) / PXDEC_Y)


def row_of_y_flat(y):
    return ROW0 - PXDEC_Y * np.log10(np.asarray(y, float))


def printed_line(re, alpha=1.1, beta=0.6):
    """The correlation drawn on the figure, (Sh-2)/Sc^(1/3) = 1.1 Re^0.6."""
    return alpha * np.asarray(re, float) ** beta


def curve_mask(shape, halfwidth):
    """Pixels within `halfwidth` of the COMPUTED position of 1.1 Re^0.6.

    `halfwidth = 0` returns an empty mask - the band is switched off - rather
    than the one-pixel line, so that the band sweep starts from "no band at
    all" and the 61 px straight-ink removal is left doing the work alone.
    """
    m = np.zeros(shape, bool)
    if halfwidth <= 0:
        return m
    cols = np.arange(int(FRAME["left"]), int(FRAME["right"]) + 1)
    rows = row_of_y(printed_line(re_of_col(cols)), cols)
    good = np.isfinite(rows) & (rows > FRAME["top"]) & (rows < FRAME["bottom"])
    m[np.round(rows[good]).astype(int), cols[good]] = True
    return ndi.binary_dilation(m, ndi.generate_binary_structure(2, 2),
                               iterations=halfwidth)


def tick_mask(ink, reach=28):
    h, w = ink.shape
    b = ink > 0.5
    lab, _ = ndi.label(b, ndi.generate_binary_structure(2, 2))
    seed = np.zeros((h, w), bool)
    for r in (int(FRAME["top"]), int(FRAME["bottom"])):
        seed[r - 3:r + 4, int(FRAME["left"]):int(FRAME["right"])] = True
    for c in (int(FRAME["left"]), int(FRAME["right"])):
        seed[int(FRAME["top"]):int(FRAME["bottom"]), c - 3:c + 4] = True
    ids = np.unique(lab[seed & b])
    touch = np.isin(lab, ids[ids > 0])
    zone = np.zeros((h, w), bool)
    zone[:int(FRAME["top"]) + reach, :] = True
    zone[int(FRAME["bottom"]) - reach:, :] = True
    zone[:, :int(FRAME["left"]) + reach] = True
    zone[:, int(FRAME["right"]) - reach:] = True
    return touch & zone


def trace_line(ink, long_len=61, n_ang=18, window=22, max_width=6):
    """Trace the ink of the drawn correlation, column by column.

    This is the figure's own control object: the curve the paper printed inside
    its own plot is known to be exactly (Sh-2)/Sc^(1/3) = 1.1 Re^0.6, so fitting
    a power law to its ink and comparing with 1.1 / 0.6 measures the axis
    calibration against something whose answer is known in advance.

    A column is used only if the straight ink near the computed position forms
    exactly ONE run no wider than `max_width` px.  Two runs, or a fat one, means
    a marker is sitting on the line there and the centroid would be pulled off
    it; that filter is what separates the line from the marker chains along it.
    Returns (columns, rows, widths).
    """
    import glyphfit as G
    straight = G.line_ink(ink, long_len, n_ang).astype(bool) & (ink > 0.5) \
        & ~tick_mask(ink)
    cols = np.arange(int(FRAME["left"]) + 30, int(FRAME["right"]) - 30)
    comp = row_of_y(printed_line(re_of_col(cols)), cols)
    oc, orow, ow = [], [], []
    for c, rc in zip(cols, comp):
        if not np.isfinite(rc) or rc < FRAME["top"] + 25 or rc > FRAME["bottom"] - 25:
            continue
        lo = int(rc) - window
        idx = np.where(straight[lo:lo + 2 * window + 1, c])[0]
        if len(idx) == 0:
            continue
        groups, cur = [], [idx[0]]
        for i in idx[1:]:
            if i - cur[-1] <= 1:
                cur.append(i)
            else:
                groups.append(cur); cur = [i]
        groups.append(cur)
        if len(groups) != 1 or len(groups[0]) > max_width:
            continue
        oc.append(float(c)); orow.append(lo + float(np.mean(groups[0])))
        ow.append(len(groups[0]))
    return np.array(oc), np.array(orow), np.array(ow, int)


SHAPES = ["circle", "square", "diamond", "tri_up", "tri_down", "tri_left",
          "tri_right", "star", "wide_dia", "hexag"]
SIZES = [13, 15, 17, 19, 21, 23]
STROKE, R, SEARCH, GUARD = 0.26, 14, 5, 0.25


def prepare(long_len=45, n_ang=18, line_hw=7):
    ink = load()
    calibrate(ink)
    h, w = ink.shape
    care = np.ones((h, w), float)
    care[:int(FRAME["top"]) + 4, :] = 0.0
    care[int(FRAME["bottom"]) - 4:, :] = 0.0
    care[:, :int(FRAME["left"]) + 4] = 0.0
    care[:, int(FRAME["right"]) - 4:] = 0.0

    import glyphfit as G
    b = ink > 0.5
    removed = (G.line_ink(ink, long_len, n_ang).astype(bool)
               | tick_mask(ink)
               | (G.line_ink(ink, 25, n_ang).astype(bool) & curve_mask((h, w), line_hw))
               ) & b
    care = care * (1.0 - removed)

    allow = np.zeros((h, w), bool)
    allow[int(FRAME["top"]) + 10:int(FRAME["bottom"]) - 10,
          int(FRAME["left"]) + 10:int(FRAME["right"]) - 10] = True
    for r0, r1, c0, c1 in TEXT_BOXES:
        allow[r0:r1, c0:c1] = False
    return b.astype(float), care, allow, removed
