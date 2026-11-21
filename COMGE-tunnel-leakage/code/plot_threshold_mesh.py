"""Reproduce the near-threshold mesh comparison in response Fig. R2.2.

The upper row contains the original far-field comparison. The lower row shows
the J7 meshes at their physical dimensions and the peak margins for the two
realizations nearest the threshold, one below it and one above it.
"""

import argparse
import csv
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.ticker import MaxNLocator
import numpy as np
from PIL import Image

BLUE = "#186585"
OCHRE = "#BC751E"
INK = "#233741"
INTERFACE = "#9c6054"


def configure_plot(font_dir):
    if font_dir:
        for path in sorted(font_dir.iterdir()):
            if path.suffix.lower() == ".ttf" and "euclid" in path.name.lower():
                font_manager.fontManager.addfont(str(path))
    try:
        font_manager.findfont("Euclid", fallback_to_default=False)
    except ValueError as error:
        raise ValueError(
            "Install the Euclid fonts or pass their directory with --fonts."
        ) from error
    plt.rcParams.update(
        {
            "font.family": "Euclid",
            "font.size": 8.5,
            "axes.linewidth": 0.64,
            "axes.unicode_minus": False,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "axes.labelpad": 2,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def draw_mesh(axis, geometry, mesh, title, spacing_label):
    axis.add_collection(
        PolyCollection(
            geometry[f"{mesh}_qTunnelGlobal"],
            facecolors="#CBD2D4",
            edgecolors="#657982",
            linewidths=0.43,
        )
    )
    axis.add_collection(
        PolyCollection(
            geometry[f"{mesh}_tSoil"],
            facecolors="#F4EEE0",
            edgecolors="#8D8269",
            linewidths=0.41,
        )
    )
    axis.add_collection(LineCollection(geometry[f"{mesh}_band"], colors=INTERFACE, linewidths=1.9))
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(-52, 52)
    axis.set_ylim(-22, 66)
    axis.set_xticks([-40, 0, 40])
    axis.set_yticks([0, 30, 60])
    axis.set_title(title, loc="left", fontsize=8.7, pad=7)
    axis.text(
        0.97,
        0.95,
        spacing_label,
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color=INK,
        bbox=dict(boxstyle="square,pad=.17", facecolor="white", edgecolor="none", alpha=0.92),
    )
    axis.set_xlabel("Tangential offset (mm)", fontsize=7.8)
    axis.tick_params(labelsize=7.2)


def draw_peak_margins(axis, rows):
    axis.set_title("(f) Peak relative to threshold", loc="left", fontsize=8.7, pad=7)
    all_margins = [float(row["margin_pct"]) for row in rows]
    bound = max(0.65, max(abs(value) for value in all_margins) * 1.2)
    axis.set_xlim(-bound, bound)
    axis.set_ylim(-0.46, 1.46)
    axis.axvline(0, color=INTERFACE, ls="--", lw=1)
    # The fixed case IDs identify the nearest peaks below and above the threshold.
    for y, case_id in [(1, "1665"), (0, "364")]:
        case_rows = [row for row in rows if row["case_id"] == case_id]
        if len(case_rows) != 2 or {row["mesh"] for row in case_rows} != {"original", "refined"}:
            raise ValueError(f"Expected one original and one refined peak for case {case_id}")
        margins = {row["mesh"]: float(row["margin_pct"]) for row in case_rows}
        original = margins["original"]
        refined = margins["refined"]
        axis.plot([original, refined], [y + 0.12, y - 0.12], c="#8b979e", lw=0.85, zorder=1)
        axis.plot(
            original, y + 0.12, "o", mec=BLUE, mfc=BLUE, ms=4.3, label="10 mm" if y == 1 else None
        )
        axis.plot(
            refined,
            y - 0.12,
            "D",
            mec=OCHRE,
            mfc="white",
            mew=1,
            ms=4.2,
            label="5 mm" if y == 1 else None,
        )
    axis.set_yticks([0, 1], ["Upper", "Lower"])
    axis.tick_params(axis="y", length=0, pad=3)
    axis.tick_params(labelsize=7.5)
    axis.xaxis.set_major_locator(MaxNLocator(4, symmetric=True))
    axis.set_xlabel("Peak margin (%)", fontsize=8)


def combine_rows(upper, lower):
    if upper.width != lower.width:
        raise ValueError("The two figure rows have different widths")
    # Trim only the verified white margins between the two 600-dpi rows.
    white_rgb = ((255, 255),) * 3
    top_margin = 120
    bottom_margin = 115
    if (
        upper.crop((0, upper.height - bottom_margin, upper.width, upper.height)).getextrema()
        != white_rgb
    ):
        raise ValueError("The upper row has marks in its bottom margin")
    if lower.crop((0, 0, lower.width, top_margin)).getextrema() != white_rgb:
        raise ValueError("The lower row has marks in its top margin")
    upper = upper.crop((0, 0, upper.width, upper.height - bottom_margin))
    lower = lower.crop((0, top_margin, lower.width, lower.height))
    figure = Image.new("RGB", (upper.width, upper.height + lower.height), "white")
    figure.paste(upper, (0, 0))
    figure.paste(lower, (0, upper.height))
    return figure


def plot(data_dir, output, font_dir=None):
    configure_plot(font_dir)
    with (data_dir / "peaks.csv").open(encoding="utf8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    fig = plt.figure(figsize=(7.16, 2.45))
    original_axis = fig.add_axes([0.076, 0.245, 0.23, 0.55])
    refined_axis = fig.add_axes([0.385, 0.245, 0.23, 0.55])
    margin_axis = fig.add_axes([0.742, 0.245, 0.235, 0.55])
    with np.load(data_dir / "local_mesh_view.npz") as geometry:
        draw_mesh(original_axis, geometry, "original", "(d) Original J7 mesh", "10 mm")
        original_axis.set_ylabel("Radial offset (mm)", fontsize=7.8)
        draw_mesh(refined_axis, geometry, "refined", "(e) Refined J7 mesh", "5 mm")
    draw_peak_margins(margin_axis, rows)
    handles, labels = margin_axis.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.855, 0.006),
        ncol=2,
        fontsize=7.6,
        frameon=False,
        handlelength=0.85,
        columnspacing=1,
        borderaxespad=0,
    )
    fig.canvas.draw()
    with io.BytesIO() as buffer:
        fig.savefig(buffer, format="png", dpi=600, facecolor="white")
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as image:
            lower = image.convert("RGB")
    with Image.open(data_dir / "far_field_panels.png") as image:
        upper = image.convert("RGB")
    output.parent.mkdir(parents=True, exist_ok=True)
    combine_rows(upper, lower).save(output, dpi=(600, 600))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data",
        type=Path,
        required=True,
        help="Directory containing local_mesh_view.npz and peaks.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Destination for the combined figure, e.g. Fig_R2_2.png",
    )
    parser.add_argument("--fonts", type=Path, help="Directory containing Euclid TrueType fonts")
    args = parser.parse_args()
    plot(args.data, args.output, args.fonts)
