"""Reproduce Fig. 12 (response Fig. R2.5) using the fitted Kriging model.

One saved Monte Carlo draw provides the hydraulic inputs for panels (a–c).
Panel (d) uses all the saved peak predictions at the same parameter setting.
"""

from dataclasses import dataclass
import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch
from scipy.spatial.distance import cdist


def configure_fonts(font_dir):
    if font_dir:
        for path in sorted(font_dir.iterdir()):
            if path.suffix.lower() != ".ttf" or "euclid" not in path.name.lower():
                continue
            font_manager.fontManager.addfont(str(path))
    try:
        for weight in ("normal", "bold"):
            font_manager.findfont(
                font_manager.FontProperties(family="Euclid", weight=weight),
                fallback_to_default=False,
            )
    except ValueError as error:
        raise ValueError(
            "Install the Euclid fonts or pass their directory with --font-dir."
        ) from error
    plt.rcParams.update(
        {
            "font.family": "Euclid",
            "font.size": 8.6,
            "axes.labelsize": 8.6,
            "xtick.labelsize": 7.7,
            "ytick.labelsize": 7.7,
            "axes.linewidth": 0.58,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Euclid",
            "mathtext.it": "Euclid:italic",
            "mathtext.bf": "Euclid:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


@dataclass
class WorkflowData:
    """Quantities used by the four panels, in the order followed by Fig. 12."""

    soil_coordinates: np.ndarray
    joint_distance: np.ndarray
    soil_log_multiplier: np.ndarray
    radius_m: np.ndarray
    training_inputs: np.ndarray
    training_peaks: np.ndarray
    correlations: np.ndarray
    predicted_peak: float
    peaks: np.ndarray
    threshold: float
    exceedances: int


def load_workflow(archive):
    """Rebuild the features and noise-free prediction for draw 1, then check the saved peak counts."""
    sys.path.insert(0, str(archive / "code/original"))
    from field_sampler import physical_fields
    from uq_workflow import physical_features

    data = archive / "data"
    inputs = (
        "data/example/geometry.npz",
        "data/figure12/PC01739_completed.json",
        "code/original/uq_workflow.py",
        "data/model/training608.npz",
        "data/model/frozen_predictor_arrays.npz",
        "data/figure12/PC01739.npz",
        "code/original/field_sampler.py",
    )
    checksums = json.loads((archive / "manifest.json").read_text(encoding="utf8"))
    for name in inputs:
        if hashlib.sha256((archive / name).read_bytes()).hexdigest() != checksums[name]:
            raise ValueError(f"Input checksum mismatch for {name}")

    with np.load(data / "example/geometry.npz") as geometry:
        soil_coordinates = geometry["soil_feature_coordinates_m"]
        joint_distance = geometry["joint_coordinates"][:, 0]
    settings = json.loads((data / "figure12/PC01739_completed.json").read_text(encoding="utf8"))[
        "parameters"
    ]

    with np.load(data / "figure12/PC01739.npz") as case:
        soil_draw = case["soil_gaussian"][0]
        joint_draw = case["joint_gaussian"][0]
        soil_log_multiplier, radius_m, _ = physical_fields(
            soil_draw, joint_draw, settings["soil_cov"], settings["joint_cov"]
        )
        geometry = {
            "nodes": soil_coordinates,
            "residual_witness_ids": np.arange(len(soil_coordinates)),
        }
        features = physical_features(
            soil_draw, joint_draw, settings["soil_cov"], settings["joint_cov"], geometry, settings
        )
        np.testing.assert_allclose(features, case["X28"][0], rtol=0.0, atol=1e-12)

        # Index 20 is the adopted threshold (eta = 1) on the 41-point grid.
        threshold_index = 20
        peaks = case["predictions"].ravel()
        threshold = float(case["thresholds"][threshold_index])
        exceedances = int((peaks > threshold).sum())
        if exceedances != int(case["exceedance_counts"][threshold_index]):
            raise ValueError("Exceedance count differs from the saved result")

    with np.load(data / "model/training608.npz") as train, np.load(
        data / "model/frozen_predictor_arrays.npz"
    ) as model:
        training_inputs = (train["X28"] - model["input_mean"]) / model["input_scale"]
        training_peaks = train["y"].reshape(-1)
        np.testing.assert_array_equal(training_inputs, model["X_train_scaled"])
        draw_inputs = (features - model["input_mean"]) / model["input_scale"]
        if np.any(np.all(draw_inputs == training_inputs, axis=1)):
            raise ValueError("The illustrated draw coincides with a training input")

        # Scale each of the 28 input distances by its fitted kernel length scale.
        length_scales = model["length_scale"]
        squared_distances = cdist(
            draw_inputs[None, :] / length_scales, training_inputs / length_scales, "sqeuclidean"
        )
        correlations = np.exp(-0.5 * squared_distances)[0]
        kernel_vector = float(model["signal_variance"]) * correlations
        # Back-transform the fitted standardized log peak into leakage units.
        predicted_peak = float(
            np.exp(model["target_mean"] + model["target_scale"] * (kernel_vector @ model["alpha"]))
        )
        np.testing.assert_allclose(predicted_peak, peaks[0], rtol=0.0, atol=1e-10)
    if (
        training_inputs.shape != (608, 28)
        or training_peaks.shape != (608,)
        or peaks.size != 131072
        or joint_distance.size != 64
    ):
        raise ValueError("Unexpected dimensions in the Fig. 12 input data")
    if threshold != 0.17827224170071532:
        raise ValueError(f"Unexpected threshold in Fig. 12 data: {threshold}")
    return WorkflowData(
        soil_coordinates,
        joint_distance,
        soil_log_multiplier,
        radius_m,
        training_inputs,
        training_peaks,
        correlations,
        predicted_peak,
        peaks,
        threshold,
        exceedances,
    )


def plot(archive, output, font_dir, alignment_tools=None):
    configure_fonts(font_dir)
    output.mkdir(parents=True, exist_ok=True)
    data = load_workflow(archive)
    density, edges = np.histogram(data.peaks, bins=55, density=True)

    ink, arrow, rule = "#173747", "#24495B", "#536773"
    water, red = "#317A93", "#AE4A35"
    earth = LinearSegmentedColormap.from_list(
        "earth_water", ["#234D72", "#3D969A", "#E8E8CB", "#C99650", "#994D35"]
    )
    fig = plt.figure(figsize=(7.16, 4.90), facecolor="white")
    left, right, width, top, bottom, height = 0.086, 0.605, 0.306, 0.607, 0.120, 0.292
    field_panel, training_panel, kernel_panel, probability_panel = [
        fig.add_axes([panel_x, panel_y, width, height])
        for panel_x, panel_y in [(left, top), (right, top), (right, bottom), (left, bottom)]
    ]
    for axis in (field_panel, training_panel, kernel_panel, probability_panel):
        axis.tick_params(direction="out", width=0.58, length=2.8, pad=2.2, colors=ink)
        for spine in axis.spines.values():
            spine.set_color(rule)
        axis.xaxis.label.set_color(ink)
        axis.yaxis.label.set_color(ink)

    def title(xpos, ypos, letter, words):
        fig.text(
            xpos, ypos, "(" + letter + ")", fontsize=10, fontweight="bold", va="baseline", color=ink
        )
        fig.text(xpos + 0.043, ypos, words, fontsize=10, va="baseline", color=ink)

    title(left, 0.947, "a", "One Monte Carlo realization")
    title(right, 0.947, "b", "Kriging training data")
    fig.text(right, 0.912, r"THM peak ($\mathrm{L\,m^{-1}\,d^{-1}}$)", fontsize=7.2, color=ink)
    title(right, 0.456, "c", "Kriging prediction")
    title(left, 0.456, "d", "Leakage exceedance probability")

    def colorbar(map_image, xpos, ypos, bar_height, ticks, labels):
        axis = fig.add_axes([xpos, ypos, 0.009, bar_height])
        bar = fig.colorbar(map_image, cax=axis, ticks=ticks)
        bar.ax.set_yticklabels(labels)
        bar.ax.tick_params(
            color=rule, labelcolor=ink, labelsize=7.1, pad=1.6, length=2.0, width=0.45
        )
        bar.outline.set_linewidth(0.45)
        bar.outline.set_edgecolor(rule)
        return bar

    # Panel (a): one physical realization at the surrogate sampling coordinates.
    field_panel.set_axis_off()
    soil_ax = fig.add_axes([left, 0.665, 0.233, 0.227])
    map_image = soil_ax.scatter(
        data.soil_coordinates[:, 0],
        data.soil_coordinates[:, 1],
        c=data.soil_log_multiplier / np.log(10),
        s=7.5,
        lw=0,
        cmap=earth,
        vmin=-0.30,
        vmax=0.30,
        clip_on=False,
    )
    soil_ax.add_patch(Circle((30, 20), 6, fc="none", ec=ink, lw=0.55))
    soil_ax.set(xlim=(0, 60), ylim=(40, 0), aspect="equal")
    soil_ax.set(
        xlabel="Horizontal coordinate (m)",
        ylabel="Depth (m)",
        xticks=[0, 20, 40, 60],
        yticks=[0, 10, 20, 30, 40],
    )
    inset = fig.add_axes([0.374, 0.666, 0.093, 0.179])
    inset.plot(data.joint_distance, data.radius_m * 1e6, color=ink, lw=0.85)
    inset.set(xlim=(0, 0.45), xticks=[0, 0.4], yticks=[1.0, 1.4, 1.8])
    inset.set_title("J7 radius (µm)", fontsize=7.7, pad=4, color=ink)
    inset.set_xlabel("Distance (m)", fontsize=7.7, labelpad=3)
    for small in (inset, soil_ax):
        small.tick_params(labelsize=7.0, pad=1.5, length=2.3, width=0.5, colors=ink)
        small.xaxis.label.set_color(ink)
        small.yaxis.label.set_color(ink)
        for spine in small.spines.values():
            spine.set_color(rule)
    bar = fig.colorbar(
        map_image,
        cax=fig.add_axes([0.198, 0.543, 0.120, 0.009]),
        orientation="horizontal",
        ticks=[-0.3, 0, 0.3],
    )
    bar.ax.tick_params(labelsize=6.9, pad=1.5, length=1.7, width=0.45, colors=ink)
    bar.outline.set_linewidth(0.45)
    fig.text(left, 0.545, r"$\log_{10}(K/\overline{K})$", va="center", fontsize=7.5, color=ink)

    # Panel (b): ordered training predictors and their matching five-year THM peaks.
    training_panel.set_axis_off()
    targets_ax = training_panel.inset_axes([0, 0.75, 1.0, 0.245])
    targets_ax.scatter(np.arange(1, 609), data.training_peaks, color=water, s=1.8, lw=0)
    targets_ax.set_xlim(0.5, 608.5)
    targets_ax.set_ylim(0, max(data.training_peaks) * 1.04)
    targets_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    targets_ax.tick_params(axis="y", labelsize=6.5, length=2.0, pad=1.8, colors=ink)
    targets_ax.spines["bottom"].set_color(rule)
    targets_ax.spines["left"].set_color(rule)
    inputs_ax = training_panel.inset_axes([0, 0, 1.0, 0.66])
    map_image = inputs_ax.imshow(
        data.training_inputs.T,
        extent=(0.5, 608.5, 28.5, 0.5),
        origin="upper",
        aspect="auto",
        cmap="RdBu_r",
        vmin=-3,
        vmax=3,
        interpolation="none",
    )
    for boundary in (7.5, 17.5):
        inputs_ax.axhline(boundary, color="white", lw=0.56)
    inputs_ax.tick_params(colors=ink, length=2.8, pad=2.2)
    inputs_ax.set(
        xlabel="Training realization",
        ylabel="Input component",
        xticks=[1, 200, 400, 608],
        yticks=[1, 7, 17, 28],
    )
    colorbar(map_image, right + width + 0.008, top, height * 0.66, [-3, 0, 3], ["−3", "0", "3"])

    # Panel (c): kernel correlation of the new input with each training input.
    kernel_panel.scatter(
        np.arange(1, len(data.correlations) + 1),
        data.correlations,
        s=3.8,
        linewidths=0,
        color=water,
        alpha=0.66,
    )
    kernel_panel.set(
        xlabel="Training realization",
        ylabel="Kernel correlation",
        xticks=[1, 200, 400, 608],
        xlim=(-10, 620),
        ylim=(0.70, 1.02),
        yticks=[0.7, 0.8, 0.9, 1.0],
    )
    kernel_panel.text(
        0.025,
        0.085,
        f"Predicted peak = {data.predicted_peak:.4f}" + r" $\mathrm{L\,m^{-1}\,d^{-1}}$",
        transform=kernel_panel.transAxes,
        fontsize=8.1,
        color=ink,
    )

    # Panel (d): the leakage distribution and its probability above the adopted threshold.
    probability_panel.stairs(density, edges, fill=True, color=water, alpha=0.74, lw=0.55)
    first_tail = np.searchsorted(edges, data.threshold, side="right") - 1
    probability_panel.stairs(
        density[first_tail:],
        np.r_[data.threshold, edges[first_tail + 1 :]],
        fill=True,
        color="#B96F50",
        lw=0,
    )
    probability_panel.axvline(data.threshold, color=red, linestyle=(0, (3.8, 2.4)), lw=1.2)
    probability_panel.text(
        data.threshold + 0.008,
        11.55,
        "Assessment\nthreshold",
        color=red,
        fontsize=7.7,
        va="top",
        linespacing=1.04,
    )
    probability_panel.text(
        0.975,
        0.55,
        f"{data.exceedances:,} of {len(data.peaks):,}\nexceed ({100*data.exceedances/len(data.peaks):.2f}%)",
        transform=probability_panel.transAxes,
        ha="right",
        color=ink,
        fontsize=7.4,
        linespacing=1.12,
    )
    probability_panel.set(
        xlim=(edges[0] - 0.005, edges[-1] + 0.007),
        ylim=(0, 13),
        xticks=[0.1, 0.2, 0.3, 0.4],
        yticks=[0, 4, 8, 12],
        xlabel=r"Five-year rolling peak ($\mathrm{L\,m^{-1}\,d^{-1}}$)",
        ylabel=r"Probability density ($\mathrm{m\,d\,L^{-1}}$)",
    )

    # A clockwise route keeps the three straight arrows in dedicated gutters.
    for start, end in [
        ((0.483, 0.752), (0.540, 0.752)),
        ((0.907, 0.550), (0.907, 0.427)),
        ((0.535, 0.266), (0.407, 0.266)),
    ]:
        fig.add_artist(
            FancyArrowPatch(
                start,
                end,
                transform=fig.transFigure,
                arrowstyle="-|>,head_length=5.0,head_width=2.45",
                mutation_scale=1.50,
                lw=1.55,
                color=arrow,
                shrinkA=0,
                shrinkB=0,
                capstyle="butt",
                joinstyle="miter",
            )
        )
    fig.canvas.draw()
    if alignment_tools:
        sys.path.insert(0, str(alignment_tools))
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig,
            axes=[field_panel, training_panel, kernel_panel, probability_panel],
            panel_ids=list("abcd"),
            row_groups=[["a", "b"], ["d", "c"]],
            column_groups=[["a", "d"], ["b", "c"]],
            json_out=output / "alignment.json",
            tolerance_pt=1.0,
        )
    fig.savefig(output / "Fig12_Kriging_workflow.png", dpi=600)
    fig.savefig(
        output / "Fig12_Kriging_workflow.pdf",
        dpi=600,
        metadata={
            "Title": "Data flow through the Kriging model",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    fig.savefig(output / "Fig12_preview.png", dpi=170)
    plt.close(fig)
    return data.training_inputs, data.correlations, density, edges


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default=Path(__file__).resolve().parents[1], type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--font-dir", type=Path, help="Directory containing Euclid TrueType fonts")
    parser.add_argument(
        "--alignment-tools",
        type=Path,
        help="Optional directory containing audit_panel_alignment.py for a panel alignment report",
    )
    args = parser.parse_args()
    plot(args.archive, args.output, args.font_dir, args.alignment_tools)
