"""Drawing the TAST plot (same style as the original notebook)."""
from __future__ import annotations

from pathlib import Path
from typing import Callable, List

import matplotlib
matplotlib.use("Agg")  #draw to files only, no window needed
import matplotlib.patches as patches  #noqa: E402
import matplotlib.pyplot as plt  #noqa: E402
import pandas as pd  #noqa: E402

from .pipeline import Row  #noqa: E402

LINE_COLOURS = ["b", "g", "m", "c", "w", "grey", "k", "silver", "beige", "teal",
                "lightskyblue", "r", "y", "lightgreen", "peru"]
REVERSE_COLOUR = "lightsteelblue"
MAX_PIXELS = 65000  #matplotlib cannot write PNGs taller or wider than 2^16 pixels


def figure_size(seq_len: int, n_rows: int):
    return min(int(seq_len / 1100), 50), max(int(n_rows / 2.5), 2)


def draw(rows: List[Row], seq_len: int, row_numbers: bool = True, dpi: int = 140):
    """Return a matplotlib figure with one panel per row.

    Red: alignment regions. Lines: extent of each transcript (one colour per transcript in
    a row). Gold outline: region overlaps an ORF. Gold fill: region overlaps the largest ORF.
    Numbers: region numbers (see layout.region_labels). Background: colinear blocks of regions,
    pink where numbers increase along the genome (or a single region), blue where they decrease.
    """
    w, h = figure_size(seq_len, len(rows))
    nrows = len(rows)
    with plt.rc_context({"agg.path.chunksize": 10000, "figure.figsize": [w, h], "figure.dpi": dpi}):
        fig, axes = plt.subplots(nrows=nrows, ncols=1, sharex=False, squeeze=False)
        axes = axes[:, 0]
        plt.subplots_adjust(top=0.9, bottom=0.1, hspace=0)
        axes[nrows - 1].set_xlabel("")

        for idx, row in enumerate(rows):
            ax = axes[idx]
            #alignment regions as red blocks
            starts = [r[0] for r in row.regions]
            ends = [r[1] for r in row.regions]
            area = pd.DataFrame()
            area["coords"] = [0] + [c for r in row.regions for c in r] + [s - 1 for s in starts] \
                + [e + 1 for e in ends] + [seq_len]
            area.sort_values(by="coords", inplace=True)
            area["coverage"] = [0] + [0, 1, 1, 0] * len(row.regions) + [0]
            area.plot.area(x="coords", y="coverage", figsize=(w, h), grid=True, legend=False, ax=ax,
                           color="r", linewidth=1, xlim=[0, seq_len], zorder=1)
            #transcript extents
            lines = pd.DataFrame()
            for i, line in enumerate(row.ranges):
                lines["coord"] = line
                lines["midline"] = 0.5
                lines.plot.line(x="coord", y="midline", figsize=(w, h), grid=True, legend=False, ax=ax,
                                color=LINE_COLOURS[i % len(LINE_COLOURS)], linewidth=4, xlim=[0, seq_len])
            #regions overlapping any ORF: gold outline
            for c, flag in zip(row.regions, row.in_orf):
                if flag:
                    ax.add_patch(patches.Rectangle((c[0], 0), c[1] - c[0], 1, color="gold", fill=None, zorder=1))
            #regions overlapping the largest ORF: gold fill
            for c, flag in zip(row.regions, row.in_largest_orf):
                if flag:
                    ax.add_patch(patches.Rectangle((c[0], 0), c[1] - c[0], 1, color="gold", zorder=99))
                    ax.add_patch(patches.Rectangle((c[0], 0), c[1] - c[0], 1, color="gold", fill=None, zorder=1))
            #region numbers
            for x, label in zip(row.label_x, row.labels):
                fig.text(x, int((nrows - idx) / (nrows + 1)) + 0.5, f"{label}", horizontalalignment="center",
                         transform=ax.transData, verticalalignment="center", fontsize=15, color="black")
            #colinear blocks of regions: pink background (forward or single region), blue (reverse)
            for c, direction in zip(row.frames, row.frame_directions):
                colour = REVERSE_COLOUR if direction == -1 else "pink"
                ax.add_patch(patches.Rectangle((c[0], 0), c[1] - c[0], 1, color=colour, zorder=0))

        def style_x(ax, **sides):
            ax.tick_params(axis="y", which="major", labelsize=30, pad=0, width=0, size=0)
            ax.tick_params(axis="y", which="minor", labelsize=0, pad=0, width=0, size=0)
            ax.tick_params(axis="x", which="major", labelsize=30, pad=15, width=3, size=8, **sides)
            ax.tick_params(axis="x", which="minor", labelsize=30, pad=15, width=2, size=5, **sides)
            ax.grid(which="major", axis="x", color="black", linestyle="--", dashes=(5, 5), linewidth=1.35, alpha=0.65)
            ax.grid(which="minor", axis="x", color="black", linestyle="-", dashes=(4, 5), linewidth=0.85, alpha=0.45)
            for side in ["top", "bottom", "left", "right"]:
                ax.spines[side].set_linewidth(2.5)
            axes[0].xaxis.tick_top()

        for idx in range(nrows):
            ax = axes[idx]
            if idx != 0 and idx != nrows - 1:  #middle rows: no x labels
                ax.set_xticklabels([])
                ax.set_xlabel("")
                ax.minorticks_on()
                style_x(ax, labelbottom=True, labeltop=True, bottom=True, top=True)
            elif idx == 0:  #top row: x labels on top
                ax.set_xlabel("")
                ax.minorticks_on()
                style_x(ax, labelbottom=0, labeltop=True, bottom=0, top=True)
            elif idx == nrows - 1:  #bottom row: x labels below
                style_x(ax, labelbottom=True, labeltop=0, bottom=True, top=0)
                ax.set_xlabel("")
                ax.minorticks_on()
            for label in ax.get_yticklabels():
                label.set_visible(False)
            ax.set_yticklabels([])
            if row_numbers:
                ax.text(-0.004, 0.5, str(idx), transform=ax.transAxes, ha="right", va="center", fontsize=22)
    return fig


PAD_INCHES = 0.3  #white border around the saved picture


def save(fig, path: Path, dpi: int, log: Callable[[str], None] = print) -> None:
    """Save the figure so that everything drawn (including the axis numbers above the top row
    and below the bottom row) fits in the picture, with empty margins trimmed. Lowers the
    resolution if the image would be too large for PNG."""
    bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    w, h = bbox.width + 2 * PAD_INCHES, bbox.height + 2 * PAD_INCHES
    if Path(path).suffix.lower() == ".png" and max(w, h) * dpi > MAX_PIXELS:
        new_dpi = int(MAX_PIXELS / max(w, h))
        log(f"Plot is too large for a {dpi} dpi PNG; saving at {new_dpi} dpi instead "
            "(use --format pdf or svg for full resolution)")
        dpi = new_dpi
    fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=PAD_INCHES)
    plt.close(fig)
