"""fairness_euler.py — Euler diagram of the fairness tiers at n = 7.

Draws how the four fairness notions nest and overlap, all sitting inside a
single outer rectangle of iso classes:
regular ⊂ balanced ⊂ inclusive (a strict nested chain), with completely-mixed
drawn as a crossing rounded rectangle (overlaps regular and balanced, spills
into the inclusive-but-not-balanced area). Not area-proportional — every region
is labelled with its exact count. Saves to plots/fairness_euler.png.

    uv run fairness_euler.py

The layout is computed: the five boxes are the only geometric constants; every
label, tier tag and region count is placed by deriving it from box edges and a
few named gaps, with a point→data-unit factor (PT) so text blocks are centred
rather than hand-nudged.

Region counts are the known n = 7 census values (the six inclusive regions sum
to inclusive = 10,525); completely-mixed exists only at odd n (parity theorem).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")

# ---- palette (restrained: blue sequential for nested tiers, amber accent for CM)
C_ISO_FILL, C_ISO_EDGE = "#f4f6f9", "#9aa7b4"
C_INC_FILL, C_INC_EDGE = "#e3edf8", "#5b8fd6"
C_BAL_FILL, C_BAL_EDGE = "#b9d3f0", "#2f6fb8"
C_REG_FILL, C_REG_EDGE = "#7fb0e0", "#1b4f8a"
C_CM_FILL,  C_CM_EDGE  = "#f2a63c", "#c9781a"
INK, INK_SOFT = "#1f2933", "#52616b"
CM_INK = "#7a3f00"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "figure.dpi": 200})

# ---- figure / axes (data coords run 0..14 x 0..10) --------------------------
FIG_W, FIG_H = 12.0, 6.6
AX = [0.02, 0.02, 0.96, 0.88]
# data-y units per typographic point: dpi cancels, so this is exact for the axes.
PT = 10.0 / (FIG_H * AX[3] * 72)

fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor("white")
axE = fig.add_axes(AX)
axE.set_xlim(0, 14); axE.set_ylim(0, 10); axE.axis("off")


class Box:
    def __init__(self, x0, y0, x1, y1):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1

    left   = property(lambda s: s.x0)
    right  = property(lambda s: s.x1)
    bottom = property(lambda s: s.y0)
    top    = property(lambda s: s.y1)
    cx     = property(lambda s: (s.x0 + s.x1) / 2)
    cy     = property(lambda s: (s.y0 + s.y1) / 2)
    w      = property(lambda s: s.x1 - s.x0)
    h      = property(lambda s: s.y1 - s.y0)


# ============================================================ geometry (the only
# hand-set numbers). CM spans past `balanced` on the right and cuts through the
# left tiers. The mid-row counts sit on a fixed MID_Y (below), independent of the
# box edges, so a box can shift without dragging the numbers with it.
iso = Box(0.30, 0.70, 12.40, 8.20)
inc = Box(0.80, 1.30, 11.90, 6.90)
bal = Box(1.30, 1.70,  8.10, 6.40)
reg = Box(1.80, 2.00,  5.70, 4.60)
cm  = Box(3.65, 1.80, 11.50, 5.20)
bal.y0 = (inc.bottom + reg.bottom) / 2   # floor midway between inc & reg floors
bal.y1 = inc.top - (bal.y0 - inc.bottom)  # top margin = bottom margin → centred in inc
MID_Y = cm.cy                          # shared centre-line for the mid-row counts

TAG_GAP  = 0.34     # gap from a box's top edge down to its tier tag
TAG_FS   = 13       # tier-tag font size
CLR      = 0.22     # clearance kept between a tag and the count beneath it
NUM_FS   = 14       # every region count, one size
LABEL_FS = 10.5     # every region sub-label, one size (fits the narrow columns)

# ============================================================ draw the nested
# boxes
axE.add_patch(Rectangle((iso.x0, iso.y0), iso.w, iso.h, linewidth=1.6,
                        edgecolor=C_ISO_EDGE, facecolor=C_ISO_FILL, zorder=0))

def rounded(b, ec, fc, z, alpha=1.0):
    axE.add_patch(FancyBboxPatch((b.x0, b.y0), b.w, b.h,
                  boxstyle="round,pad=0,rounding_size=0.45", linewidth=1.8,
                  edgecolor=ec, facecolor=fc, alpha=alpha, zorder=z))

rounded(inc, C_INC_EDGE, C_INC_FILL, 1)
rounded(bal, C_BAL_EDGE, C_BAL_FILL, 2)
rounded(reg, C_REG_EDGE, C_REG_FILL, 3)
rounded(cm,  C_CM_EDGE,  C_CM_FILL,  4, alpha=0.42)   # translucent so overlaps blend
rounded(cm,  C_CM_EDGE,  "none",     6)               # crisp opaque outline on top

# ============================================================ iso header (top
# band), stacked with real point-height spacing
h_title_y = iso.top - 0.20
axE.text(0.70, h_title_y, "iso classes", ha="left", va="top",
         fontsize=15, fontweight="bold", color=INK, zorder=8)
axE.text(0.70, h_title_y - 15 * PT - 0.16,
         "2,142,288 total, only 10,525 inclusive", ha="left", va="top",
         fontsize=11.5, style="italic", color=INK_SOFT, zorder=8)

# ============================================================ tier tags — each a
# fixed gap below its own box top
def tag(text, x, box, color, ha):
    axE.text(x, box.top - TAG_GAP, text, ha=ha, va="center",
             fontsize=TAG_FS, fontweight="bold", color=color, zorder=7)

colA = (reg.left + cm.left) / 2      # regular-only
colD = (bal.right + cm.right) / 2    # completely-mixed-only
tag("inclusive",        inc.right - 0.20, inc, C_INC_EDGE, "right")
tag("balanced",         bal.left + 0.35,  bal, C_BAL_EDGE, "left")
tag("regular",          colA,             reg, C_REG_EDGE, "center")
tag("completely mixed", colD,             cm,  C_CM_EDGE,  "center")

# ============================================================ count + label
# helpers.  `count_at` puts the number's centre at num_y and hangs the label
# beneath it; `centred_in` centres the whole number+label block in a band.
LABEL_GAP = 0.35 * LABEL_FS * PT     # count-to-label gap
LINE_H    = LABEL_FS * PT * 1.08     # one label line, data units

def _draw(cx, num_y, count, label, color):
    axE.text(cx, num_y, f"{count:,}", ha="center", va="center",
             fontsize=NUM_FS, fontweight="bold", color=color, zorder=8)
    axE.text(cx, num_y - NUM_FS * PT / 2 - LABEL_GAP, label, ha="center",
             va="top", fontsize=LABEL_FS, color=color, zorder=8, linespacing=1.08)

def count_at(cx, num_y, count, label, color=INK):
    _draw(cx, num_y, count, label, color)

def centred_in(cx, y_lo, y_hi, count, label, color=INK):
    n = label.count("\n") + 1
    block = NUM_FS * PT + LABEL_GAP + n * LINE_H
    num_y = (y_lo + y_hi) / 2 + block / 2 - NUM_FS * PT / 2
    _draw(cx, num_y, count, label, color)

# --- mid row: all four counts share MID_Y (a clean aligned row) and are centred
# horizontally in their column strips
colB = (cm.left + reg.right) / 2     # regular ∩ CM
colC = (reg.right + bal.right) / 2   # balanced ∩ CM
count_at(colA, MID_Y, 3,    "regular\nonly")
count_at(colB, MID_Y, 10,   "regular &\ncompletely\nmixed")
count_at(colC, MID_Y, 92,   "balanced &\ncompletely\nmixed")
count_at(colD, MID_Y, 7166, "completely mixed,\nnot balanced", color=CM_INK)

# --- top-of-box cells: centred in the band between the tier tag and the tier
# below, so they never crowd the tag
bal_band_top = (bal.top - TAG_GAP) - TAG_FS * PT / 2 - CLR
inc_band_top = (inc.top - TAG_GAP) - TAG_FS * PT / 2 - CLR
centred_in((bal.left + cm.left) / 2, reg.top, bal_band_top, 70, "balanced only")
centred_in((bal.right + inc.right) / 2, cm.top, inc_band_top, 3184,
           "inclusive only", color=C_INC_EDGE)

# ============================================================ figure title,
# centred over the iso rectangle
iso_cx_fig = AX[0] + (iso.cx / 14) * AX[2]
iso_top_fig = AX[1] + (iso.top / 10) * AX[3]
fig.text(iso_cx_fig, iso_top_fig + 0.05,
         "Fair Rock-Paper-Scissors games at n = 7",
         ha="center", va="center", fontsize=18, fontweight="bold", color=INK)

out = os.path.join(PLOTS_DIR, "fairness_euler.png")
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
print(f"saved {out}")
