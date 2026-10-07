"""Fig. 1 of the paper: plan of the L-shaped test area, anchors and the measured test points.

    python plot_fig1.py          # -> figures/fig1_test_area.png and .pdf

Test points (64 interior + 4 exterior control points) are read from the raw CSV files (ground-truth
``gt_x``, ``gt_y``). Room dimensions, the partition and the corridor are from Sec. IV of the paper; pillar
positions were measured from the published figure (about +/- 2 cm).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
from scipy.spatial import ConvexHull

from uwb import ANCHORS, ROOT, load

MAIN, SECOND, PARTITION, CORRIDOR = (5.965, 19.120), (5.461, 6.125), 0.14, 1.70
SECOND_XY = (MAIN[0] + PARTITION, 13.10)                       # top-left corner of the secondary room
PILLARS = {"S1": (1.09, 1.30), "S2": (5.50, 1.30), "S3": (1.09, 11.75), "S4": (5.50, 11.75)}
PILLAR_SIZE = (0.77, 0.87)
WALL, GREEN, RED = "#1a1a1a", "#2f6b3d", "#aa3333"
PILLAR_FACE, PILLAR_EDGE = "#d9d9d9", "#595959"


def test_points():
    """Distinct ground-truth positions of all measured test points."""
    pts = load("MK8000").drop_duplicates("point_id")[["gt_x", "gt_y"]].to_numpy()
    return pts


def draw(outdir: Path = ROOT / "figures"):
    plt.rcParams.update({"font.family": "DejaVu Sans", "hatch.linewidth": 1.7})
    fig = plt.figure(figsize=(6.97, 11.15))                    # 137 px/m at 300 dpi
    ax = fig.add_axes([0.13, 0.04, 5.92 / 6.97, 9.45 / 11.15])

    # anchor polygon (convex hull)
    xy = np.array([a[:2] for a in ANCHORS.values()])
    hull = xy[ConvexHull(xy).vertices]
    ax.add_patch(Polygon(hull, facecolor="#4682b4", alpha=0.14, edgecolor="none", zorder=1))
    ax.add_patch(Polygon(hull, fill=False, edgecolor="#3b6e9c", alpha=0.55, lw=1.6, zorder=1.5))

    # rooms, partition, corridor
    (mw, ml), (sw, sl), (sx, sy) = MAIN, SECOND, SECOND_XY
    ax.add_patch(Rectangle((0, 0), mw, ml, fill=False, edgecolor=WALL, lw=3.5, zorder=4))
    ax.add_patch(Rectangle((mw, sy), PARTITION, sl, facecolor=WALL, edgecolor=WALL, lw=3.5, zorder=4))
    ax.add_patch(Rectangle((sx, sy), sw, sl, fill=False, edgecolor=WALL, lw=3.5, zorder=4))
    ax.add_patch(Rectangle((sx, sy - CORRIDOR), sw, CORRIDOR, facecolor="#f0f0f0", edgecolor="#8c8c8c",
                           lw=2.2, linestyle=(0, (3, 2.2)), zorder=2))
    ax.text(sx + sw / 2, sy - CORRIDOR / 2, "corridor 1.70 m", ha="center", va="center", fontsize=12,
            color="#505050", zorder=3)
    ax.text(sx + 0.35, 16.65, "14 cm partition", rotation=90, ha="center", va="center", fontsize=10.8, zorder=6)

    # structural pillars
    pw, ph = PILLAR_SIZE
    for name, (cx, cy) in PILLARS.items():
        ax.add_patch(Rectangle((cx - pw / 2, cy - ph / 2), pw, ph, facecolor=PILLAR_FACE, edgecolor=PILLAR_EDGE,
                               hatch="/", lw=1.7, zorder=3))
        ax.text(cx, cy, name, ha="center", va="center", fontsize=11, color="#333333", zorder=6)

    # measured test points, split by the anchor hull
    pts = test_points()
    inside = np.array([MplPath(hull).contains_point(p) for p in pts])
    ax.plot(*pts[inside].T, "o", ms=3.3, color=GREEN, mew=0, ls="none", zorder=5)
    if (~inside).any():
        ax.plot(*pts[~inside].T, "o", ms=3.1, mfc="white", mec=RED, mew=1.1, ls="none", zorder=5)

    # anchors and labels (left edge / baseline measured from the published figure)
    labels = {"A0": (0.58, -0.356), "A1": (4.19, -0.356), "A2": (2.05, 11.28), "A3": (3.50, 18.37),
              "A4": (8.52, 15.34)}
    for name, (x, y, z) in ANCHORS.items():
        ax.plot([x], [y], "s", ms=12.7, mfc="white", mec=WALL, mew=3.0, zorder=8)
        lx, ly = labels[name]
        ax.text(lx, ly, f"{name} ({z:.2f} m)", ha="left", va="baseline", fontsize=12, zorder=9,
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85))

    ax.set_xlim(-0.62, 12.35)
    ax.set_ylim(19.92, -0.78)
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(range(0, 13, 2))
    ax.set_yticks(range(0, 19, 2))
    ax.tick_params(labelsize=11.5, length=4.5, width=1.3, direction="out")
    ax.set_xlabel("x [m]", fontsize=14)
    ax.set_ylabel("y [m]", fontsize=14)

    handles = [Line2D([], [], marker="s", ls="none", ms=11.4, mfc="white", mec=WALL, mew=3, label="anchor"),
               Line2D([], [], marker="o", ls="none", ms=9, mfc=GREEN, mec=GREEN, label="test point, inside hull")]
    if (~inside).any():
        handles.append(Line2D([], [], marker="o", ls="none", ms=6.6, mfc="white", mec=RED, mew=2.4,
                              label="test point, outside hull"))
    handles += [Patch(facecolor="#cfe0ee", edgecolor="#3b6e9c", lw=3, label="anchor polygon"),
                Patch(facecolor=PILLAR_FACE, edgecolor=PILLAR_EDGE, lw=3, label="structural pillar")]
    leg = ax.legend(handles=handles, loc="upper right", fontsize=11.3, framealpha=1, fancybox=True,
                    edgecolor="#bebebe", borderpad=0.5, labelspacing=0.55, handlelength=1.2, handleheight=1.2,
                    handletextpad=0.9, borderaxespad=0, bbox_to_anchor=(11.6, -0.57), bbox_transform=ax.transData)
    leg.set_zorder(20)
    leg.get_frame().set_linewidth(2.5)

    outdir.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"fig1_test_area.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.1, facecolor="white")
    plt.close(fig)
    print(f"{len(pts)} test points ({inside.sum()} inside the hull) -> {(outdir / 'fig1_test_area.png').relative_to(ROOT)}")


if __name__ == "__main__":
    plt.switch_backend("Agg")                                  # headless
    draw()
