"""Regenerate both A3.4 datasets and the review overlays, deterministically."""
import json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

HERE = str(Path(__file__).parent)
sys.path.insert(0, HERE)
import setup2 as F2, setup3 as F3, glyphfit as G, overlay as OV

REPO = str(Path(__file__).parent.parent)
DATA = f"{REPO}/page/data"
REVIEW = f"{REPO}/review"

FIG2_AUDIT_REJECT = [7, 9]     # visual audit: #7 sits on the "Single spheres" leader,
                               # #9 on the eq. (11) line where the leader crosses it
FIG2_FIT_OK = 0.50

# ---------------------------------------------------------------- figure 2
print("Figure 2 ------------------------------------------------------")
ink2, care2, allow2, rem2 = F2.prepare()
prev = pd.read_csv(f"{DATA}/wakao-funazkri-1978-fig2.csv")
if "Re_prev" in prev.columns:                     # already emitted once
    prev = prev[prev["pass"] == "2026-07-30"][["Re_prev", "Sh_prev"]]
    prev.columns = ["Re", "Sh"]
ys = F2.row_of_sh(prev.Sh.values)
xs = F2.col_of_re(prev.Re.values)
t = time.time()
rev = []
for a, b in zip(ys, xs):
    g = G.fit_glyph(ink2, care2, float(a), float(b), F2.SHAPES, F2.SIZES, False,
                    F2.STROKE, F2.R, search=F2.SEARCH, guard=F2.GUARD)
    g.update(old_y=float(a), old_x=float(b), origin="reviewed-2026-07-30")
    rev.append(g)
allm, _ = G.pursue_batch(ink2, care2, allow2, F2.SHAPES, F2.SIZES, filled=False,
                         stroke=F2.STROKE, R=F2.R, thresh=0.40, keep_thresh=0.55,
                         min_sep=12.0, rounds=6, seeds=rev, search=F2.SEARCH,
                         guard=F2.GUARD, ev_min=0.42, verbose=False)
new = [g for g in allm[len(rev):]]
kept = [g for i, g in enumerate(new) if i not in FIG2_AUDIT_REJECT]
print(f"  {len(rev)} reviewed re-centred, {len(new)} new found, "
      f"{len(new)-len(kept)} rejected by audit  ({time.time()-t:.0f} s)")

rows = []
for g in rev:
    ok = g["score"] >= FIG2_FIT_OK
    rows.append(dict(Re=float(F2.re_of_col(g["x"] if ok else g["old_x"])),
                     Sh=float(F2.sh_of_row(g["y"] if ok else g["old_y"])),
                     Re_prev=float(F2.re_of_col(g["old_x"])),
                     Sh_prev=float(F2.sh_of_row(g["old_y"])),
                     col_px=round(g["x"] if ok else g["old_x"], 2),
                     row_px=round(g["y"] if ok else g["old_y"], 2),
                     col_px_prev=round(g["old_x"], 2), row_px_prev=round(g["old_y"], 2),
                     shape_fitted=g["shape"], fit_score=round(g["score"], 3),
                     recentred=bool(ok)))
    rows[-1]["pass"] = "2026-07-30"
for g in kept:
    rows.append(dict(Re=float(F2.re_of_col(g["x"])), Sh=float(F2.sh_of_row(g["y"])),
                     Re_prev=np.nan, Sh_prev=np.nan,
                     col_px=round(g["x"], 2), row_px=round(g["y"], 2),
                     col_px_prev=np.nan, row_px_prev=np.nan,
                     shape_fitted=g["shape"], fit_score=round(g["score"], 3),
                     recentred=True))
    rows[-1]["pass"] = "2026-08-02"
df2 = pd.DataFrame(rows).sort_values(["Re", "Sh"]).reset_index(drop=True)

# Fitting a shape rather than a crosshair exposes double counts: several of the
# 2026-07-30 crosshairs turn out to be two or three marks on ONE glyph, and the
# fits then land on the same point.  Group anything within DEDUP px and give the
# group one marker_id; the CSV keeps a row per original crosshair so the effect
# of the re-centring can be measured on a paired sample, and the shipped point
# set is the unique marker_ids.
DEDUP = 6.0
ids, centres = [], []
for r in df2.itertuples():
    hit = next((k for k, (cy, cx) in enumerate(centres)
                if (r.row_px - cy) ** 2 + (r.col_px - cx) ** 2 < DEDUP ** 2), None)
    if hit is None:
        centres.append((r.row_px, r.col_px)); hit = len(centres) - 1
    ids.append(hit)
df2.insert(0, "marker_id", ids)
for c in ("Re", "Sh", "Re_prev", "Sh_prev"):
    df2[c] = df2[c].round(4)
df2.to_csv(f"{DATA}/wakao-funazkri-1978-fig2.csv", index=False)
n_uniq = df2.marker_id.nunique()
n_prev_uniq = df2[df2["pass"] == "2026-07-30"].marker_id.nunique()
print(f"  wrote {len(df2)} rows / {n_uniq} distinct glyphs "
      f"({len(rev)} old crosshairs -> {n_prev_uniq} glyphs, "
      f"{len(rev)-n_prev_uniq} were double counts)")

# Where the three printed curves' INK sits relative to their computed position.
# All three the same way is an ordinate origin error, not a scatter.
cv2 = F2.trace_curves(F2.load())
print("  drawn-ink offset from the computed curve (-ve = ink above computed):")
for k, (n, m, s) in cv2.items():
    print(f"    {k:5s} n={n:4d}  {m:+.2f} +/- {s:.2f} px"
          f"   = {100*(10**(-m/F2.PXDEC_Y)-1):+.2f} % in Sh")

# ---------------------------------------------------------------- figure 3
# The erasure band around the computed position of 1.1 Re^0.6 is the one
# mechanism in this pipeline that can destroy a real marker, and it moves the
# answer by more than any statistical standard error on the page.  So it is not
# run once at a chosen setting and disclosed in prose: the whole sweep is
# extracted and SHIPPED, one block of rows per half-width, and the page computes
# the systematic uncertainty from it.  BAND_SHIPPED is the point set the page
# uses for its headline numbers.
BAND_HW = (0, 3, 5, 7, 10, 14)
BAND_SHIPPED = 5

print("Figure 3 ------------------------------------------------------")
blocks, f3 = [], None
for hw in BAND_HW:
    ink3, care3, allow3, rem3 = F3.prepare(long_len=61, line_hw=hw)
    t = time.time()
    g3, _ = G.pursue_batch(ink3, care3, allow3, F3.SHAPES, F3.SIZES, filled=True,
                           stroke=F3.STROKE, R=F3.R, thresh=0.35, keep_thresh=0.42,
                           min_sep=9.0, rounds=10, search=F3.SEARCH, guard=F3.GUARD,
                           ev_min=0.30, verbose=False)
    if hw == BAND_SHIPPED:
        f3 = g3
    blocks.append(pd.DataFrame(
        [dict(band_hw=hw,
              Re=round(float(F3.re_of_col(g["x"])), 4),
              Sh_minus_2_over_Sc13=round(float(F3.y_of_row(g["y"], g["x"])), 4),
              col_px=round(g["x"], 2), row_px=round(g["y"], 2),
              shape_fitted=g["shape"], fit_score=round(g["score"], 3))
         for g in g3]).sort_values("Re").reset_index(drop=True))
    print(f"  band half-width {hw:2d}: {len(g3):3d} markers ({time.time()-t:.0f} s)")
df3 = pd.concat(blocks, ignore_index=True)
df3.to_csv(f"{DATA}/wakao-funazkri-1978-fig3.csv", index=False)
print(f"  wrote {len(df3)} rows ({len(BAND_HW)} band settings) to "
      f"wakao-funazkri-1978-fig3.csv; shipped set = band_hw {BAND_SHIPPED}, "
      f"{int((df3.band_hw == BAND_SHIPPED).sum())} markers")

# the drawn correlation's own ink: the calibration control
lc, lr, lw = F3.trace_line(F3.load())
dfl = pd.DataFrame(dict(col_px=np.round(lc, 2), row_px=np.round(lr, 2),
                        run_px=lw))
dfl["Re"] = np.round(F3.re_of_col(lc), 4)
dfl["Sh_minus_2_over_Sc13"] = np.round(F3.y_of_row(lr, lc), 4)
dfl.to_csv(f"{DATA}/wakao-funazkri-1978-fig3-line.csv", index=False)
print(f"  wrote {len(dfl)} traced columns of the drawn line to "
      f"wakao-funazkri-1978-fig3-line.csv")

print(f"  calibration: Re=1 at col {F3.COL0:.2f}, {F3.PXDEC_X:.2f} px/decade")
print(f"    y=1 at row {F3.ROW0_L:.2f} (col {F3.TICKCOL_L:.0f}) and "
      f"{F3.ROW0_R:.2f} (col {F3.TICKCOL_R:.0f}); "
      f"{F3.PXDEC_Y_L:.2f} / {F3.PXDEC_Y_R:.2f} px per decade")
print(f"    -> ordinate skew {F3.ROW0_R - F3.ROW0_L:+.2f} px across the plot; "
      f"the flat average would be y=1 at row {F3.ROW0:.2f}, {F3.PXDEC_Y:.2f} px/dec")

# ---------------------------------------------------------------- overlays
print("review overlays ----------------------------------------------")
img2 = np.asarray(Image.open(F2.PATH).convert("L"))
glyph2, seen_ids = [], set()
for i, r in df2.iterrows():
    if r.marker_id in seen_ids: continue
    seen_ids.add(r.marker_id)
    glyph2.append(dict(y=r.row_px, x=r.col_px, shape=r.shape_fitted, size=20.0,
                       score=r.fit_score,
                       old_y=(None if np.isnan(r.row_px_prev) else r.row_px_prev),
                       old_x=(None if np.isnan(r.col_px_prev) else r.col_px_prev),
                       origin=("new" if r["pass"] == "2026-08-02" else "reviewed")))
uniq2 = df2.drop_duplicates("marker_id")
for g, src in zip(glyph2, uniq2.itertuples()):
    g["size"] = float(next((x["size"] for x in rev + kept
                            if abs(x["x"] - src.col_px) < 0.01
                            and abs(x["y"] - src.row_px) < 0.01), 20.0))
OV.full_overlay(img2, glyph2, f"{REVIEW}/A3.4-fig2-recentred-overlay.png",
                f"Figure 2, {len(glyph2)} markers: green = shape-fitted centre, "
                f"orange x = the 2026-07-30 centre, red = found on 2026-08-02",
                dpi=200)
OV.contact_sheet(img2, glyph2, f"{REVIEW}/A3.4-fig2-recentred-audit.png", ncol=10,
                 half=28, dpi=200,
                 title="A3.4 Figure 2: every marker at 2x, numbered as in the CSV; "
                       "green outline = fitted glyph, orange x = old centre")

img3 = np.asarray(Image.open(F3.PATH).convert("L"))
glyph3 = [dict(y=g["y"], x=g["x"], shape=g["shape"], size=g["size"],
               score=g["score"], origin="new") for g in f3]
OV.full_overlay(img3, glyph3, f"{REVIEW}/A3.4-fig3-extraction-overlay.png",
                f"Figure 3, {len(glyph3)} shape-fitted markers "
                f"(erasure band half-width {BAND_SHIPPED} px)", dpi=200,
                show_old=False, numbers=False)
order = sorted(range(len(glyph3)), key=lambda i: glyph3[i]["score"])
OV.contact_sheet(img3, [glyph3[i] for i in order],
                 f"{REVIEW}/A3.4-fig3-candidate-audit.png", ncol=12, half=22,
                 dpi=200, show_old=False,
                 title="A3.4 Figure 3: every candidate at 2x, worst fit first")
json.dump(dict(fig2_n=len(df2), fig3_rows=len(df3),
               fig3_n=int((df3.band_hw == BAND_SHIPPED).sum()),
               fig3_band_n={int(h): int((df3.band_hw == h).sum()) for h in BAND_HW},
               fig3_line_n=len(dfl),
               fig3_cal=dict(col0=F3.COL0, pxdec_x=F3.PXDEC_X,
                             row0_left=F3.ROW0_L, row0_right=F3.ROW0_R,
                             pxdec_y_left=F3.PXDEC_Y_L, pxdec_y_right=F3.PXDEC_Y_R,
                             tickcol_left=F3.TICKCOL_L, tickcol_right=F3.TICKCOL_R,
                             row0_flat=F3.ROW0, pxdec_y_flat=F3.PXDEC_Y)),
          open(f"{HERE}/emit_summary.json", "w"), indent=1)
print("done")
