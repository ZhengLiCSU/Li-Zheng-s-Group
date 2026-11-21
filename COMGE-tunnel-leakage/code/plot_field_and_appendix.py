"""Render manuscript Fig. 8 from the archived numerical realization, and prepare Figs. A.1 and A.2."""

from pathlib import Path
import argparse
import json
import sys
import hashlib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
import numpy as np
from scipy.interpolate import griddata
import pymupdf as pdf


def read_field_data(archive):
    """Reconstruct the published conductivity and J7 radius realization."""
    sys.path.insert(0, str(archive / "code/original"))
    from field_sampler import physical_fields

    inputs = {
        "geometry": "data/example/geometry.npz",
        "gaussian": "data/example/adaptive_candidate_000.npz",
        "settings": "data/example/adaptive_candidate_000.json",
        "transform": "code/original/field_sampler.py",
    }
    checksums = json.loads((archive / "manifest.json").read_text(encoding="utf8"))
    for name in inputs.values():
        if hashlib.sha256((archive / name).read_bytes()).hexdigest() != checksums[name]:
            raise ValueError(f"Input checksum mismatch for {name}")
    with np.load(archive / inputs["geometry"]) as geometry:
        coordinates = np.vstack((geometry["nodes"], geometry["quadrature"]))
        joint_distance = geometry["joint_coordinates"][:, 0]
    settings = json.loads((archive / inputs["settings"]).read_text(encoding="utf8"))["parameters"]
    with np.load(archive / inputs["gaussian"]) as draw:
        log_multiplier, radius_m, ensemble_radius_m = physical_fields(
            draw["soil_gaussian"],
            draw["joint_gaussian"],
            settings["soil_cov"],
            settings["joint_cov"],
        )
    if coordinates.shape != (48479, 2) or len(joint_distance) != 64 or len(radius_m) != 64:
        raise ValueError("Unexpected dimensions for the archived soil or J7 field")
    if settings["soil_theta_m"] != [24, 6] or not np.isclose(settings["joint_theta_m"], 0.30):
        raise ValueError("The field parameters do not match the Fig. 8 realization")

    # Interpolation is used for the colour map only. The published J7 series uses the sampled values.
    grid_x, grid_y = np.meshgrid(np.linspace(0, 60, 601), np.linspace(0, 40, 401))
    field = griddata(coordinates, log_multiplier / np.log(10), (grid_x, grid_y), method="linear")
    field[(grid_x - 30) ** 2 + (grid_y - 20) ** 2 < 6**2] = np.nan
    return field, joint_distance, radius_m, ensemble_radius_m


def make_field(archive, out, fonts, align=None):
    field, joint_distance, radius_m, ensemble_radius_m = read_field_data(archive)
    out.mkdir(parents=True, exist_ok=True)
    if fonts:
        for path in sorted(fonts.iterdir()):
            if path.suffix.lower() == ".ttf" and "euclid" in path.name.lower():
                font_manager.fontManager.addfont(str(path))
    try:
        font_manager.findfont("Euclid", fallback_to_default=False)
    except ValueError as error:
        raise ValueError(
            "Install the Euclid fonts or pass their directory with --font-dir."
        ) from error
    plt.rcParams.update(
        {
            "font.family": "Euclid",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Euclid",
            "mathtext.it": "Euclid:italic",
            "mathtext.bf": "Euclid:bold",
            "svg.fonttype": "none",
        }
    )

    ink = "#233B49"
    teal = "#216582"
    ochre = "#AB641F"
    cmap = LinearSegmentedColormap.from_list(
        "hydraulic", ["#234D72", "#3D969A", "#E8E8CB", "#C99650", "#994D35"]
    )
    figure = plt.figure(figsize=(7.16, 2.60))
    soil_axis, radius_axis = [
        figure.add_axes([soil_box, 0.215, 0.327, 0.60]) for soil_box in (0.078, 0.630)
    ]
    image = soil_axis.imshow(
        field,
        extent=(0, 60, 40, 0),
        origin="upper",
        aspect="equal",
        interpolation="none",
        cmap=cmap,
        vmin=-1,
        vmax=1,
    )
    soil_axis.add_patch(Circle((30, 20), 6, facecolor="#DFE5E8", edgecolor=ink, lw=0.52))
    soil_axis.add_patch(Circle((30, 20), 5.55, facecolor="white", edgecolor=ink, lw=0.52))
    soil_axis.set(
        xlabel="Horizontal coordinate (m)",
        ylabel="Depth (m)",
        xticks=[0, 15, 30, 45, 60],
        yticks=[0, 10, 20, 30, 40],
    )
    color_axis = figure.add_axes([0.427, 0.241, 0.012, 0.548])
    colorbar = figure.colorbar(image, cax=color_axis, ticks=[-1, -0.5, 0, 0.5, 1])
    colorbar.set_label(r"$\log_{10}(K/\overline{K})$", labelpad=3)
    colorbar.ax.tick_params(labelsize=7.2, pad=2, length=2.4, width=0.55)
    radius_axis.plot(
        joint_distance,
        1e6 * radius_m,
        lw=1.25,
        marker="o",
        ms=2.3,
        color=teal,
        label="J7 radius scale",
    )
    radius_axis.axhline(
        1e6 * ensemble_radius_m,
        lw=1.05,
        ls=(0, (4, 2)),
        color=ochre,
        label="Ensemble arithmetic mean",
    )
    radius_axis.set(
        xlabel="Distance along J7 pathway (m)",
        ylabel="Radius scale (µm)",
        xlim=(0, 0.45),
        xticks=[0, 0.1, 0.2, 0.3, 0.4],
        yticks=[0.8, 1.0, 1.2, 1.4],
        ylim=(0.64, 1.49),
    )
    radius_axis.legend(
        frameon=False, loc="upper left", fontsize=7, handlelength=2.2, labelspacing=0.4
    )
    for axis, title_x, letter, title, conditions in [
        (
            soil_axis,
            0.078,
            "(a)",
            "Soil hydraulic conductivity",
            "COV = 0.75, SOF = 24 m (horizontal), 6 m (vertical)",
        ),
        (radius_axis, 0.630, "(b)", "Segment-joint radius scale", "COV = 0.22, SOF = 0.30 m"),
    ]:
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(direction="out", length=2.8, width=0.6, pad=2)
        figure.text(title_x, 0.943, letter, fontsize=9.7, weight="bold", color=ink)
        figure.text(title_x + 0.040, 0.943, title, fontsize=9.7, color=ink)
        figure.text(title_x, 0.855, conditions, fontsize=7.3, color=ink)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    soil_box = soil_axis.get_window_extent(renderer)
    colorbox = color_axis.get_tightbbox(renderer)
    radius_box = radius_axis.get_tightbbox(renderer)
    assert soil_box.x1 < colorbox.x0 < colorbox.x1 < radius_box.x0, (soil_box, colorbox, radius_box)
    if align:
        sys.path.insert(0, str(align))
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            figure,
            axes=[soil_axis, radius_axis],
            panel_ids=["a", "b"],
            row_groups=[["a", "b"]],
            json_out=out / "alignment8.json",
            tolerance_pt=1.0,
        )
    figure.savefig(out / "Fig8.pdf", dpi=600, metadata={"CreationDate": None, "ModDate": None})
    figure.savefig(out / "Fig8.png", dpi=600)
    plt.close(figure)


def appendix(src, dst, marker):
    """Remove the separate production note from the supplied vector PDF; preserve the plotted curves."""
    scale_factor = 3
    scale = pdf.Matrix(scale_factor, scale_factor)
    with pdf.open(src) as document:
        page = document[0]
        original_text = page.get_text()
        footer_blocks = [
            block
            for block in page.get_text("blocks")
            if block[1] > 0.90 * page.rect.height and marker in block[4]
        ]
        if len(footer_blocks) != 1:
            raise ValueError(f"Expected one separate footer containing {marker!r} in {src.name}")
        footer = footer_blocks[0]
        box = pdf.Rect(footer[:4])
        # Restrict text removal to the footer baseline; keep the neighbouring axis labels.
        strip = pdf.Rect(
            box.x0 - 0.4, box.y0 + 0.08 * box.height, box.x1 + 0.4, box.y1 - 0.12 * box.height
        )
        before = page.get_pixmap(matrix=scale, alpha=False)
        page.add_redact_annot(strip, fill=False)
        page.apply_redactions(images=0, graphics=0, text=0)
        expected_text = original_text.replace(footer[4], "")
        if page.get_text().split() != expected_text.split() or len(page.get_drawings()) <= 20:
            raise ValueError(
                f"Appendix text or vector paths changed outside the footer: {src.name}"
            )
        document.save(dst.with_suffix(".pdf"), garbage=4, deflate=True)

    with pdf.open(dst.with_suffix(".pdf")) as result:
        page = result[0]
        after = page.get_pixmap(matrix=scale, alpha=False)
        original_pixels = np.frombuffer(before.samples, dtype=np.uint8).reshape(
            before.height, before.width, 3
        )
        revised_pixels = np.frombuffer(after.samples, dtype=np.uint8).reshape(
            after.height, after.width, 3
        )
        changed_rows, _ = np.where(np.any(original_pixels != revised_pixels, axis=2))
        if (
            not changed_rows.size
            or changed_rows.min() < int(box.y0 * scale_factor) - 2
            or changed_rows.max() > int(box.y1 * scale_factor) + 2
        ):
            raise ValueError(f"Unexpected location of raster differences in {src.name}")
        page.get_pixmap(dpi=600, alpha=False).save(dst.with_suffix(".png"))
        return [int(changed_rows.min()), int(changed_rows.max()), int(changed_rows.size)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--appendix-originals", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--font-dir", type=Path)
    parser.add_argument(
        "--alignment-tools",
        type=Path,
        help="Optional directory containing audit_panel_alignment.py for a panel alignment report",
    )
    args = parser.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    make_field(args.archive, args.output, args.font_dir, args.alignment_tools)
    original_figures = args.appendix_originals or args.archive / "data/appendix_figures"
    for tag, name, marker in [
        ("A1", "f2_appendix_probability_mean", "0.17827224"),
        ("A2", "f2_appendix_contrast_validation_v3", "Error"),
    ]:
        print(
            tag, appendix(original_figures / (name + ".pdf"), args.output / ("Fig" + tag), marker)
        )
