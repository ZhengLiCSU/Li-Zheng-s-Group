"""Recompute the J7 mesh-refinement comparison for the two peaks nearest the threshold.

For each realization and mesh, print the five-year maximum rolling 24-hour
average leakage (L m⁻¹ d⁻¹) and its exceedance classification (0 or 1).
"""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/revision/near_threshold_mesh"
HORIZON = 1825.0  # days


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def nearest_cases(pool, threshold):
    """Select the nearest peak below and above the threshold in the original 688 cases."""
    if len(pool) != 688 or len({r["sample_id"] for r in pool}) != 688:
        raise ValueError("Expected 688 distinct physical cases.")
    peaks = [(float(r["peak_L_m_d"]), r["sample_id"]) for r in pool]
    below = max((pair for pair in peaks if pair[0] < threshold), key=lambda pair: pair[0])
    above = min((pair for pair in peaks if pair[0] > threshold), key=lambda pair: pair[0])
    ids = below[1], above[1]
    if ids != ("1665", "364"):
        raise ValueError(f"The threshold-adjacent cases are {ids}; expected 1665 and 364.")
    return ids


def original_mesh_digest(folder):
    """Link the published mesh to the original binary fingerprint saved with the field."""
    mesh_name = "data/revision/near_threshold_mesh/mesh_original.mat"
    record = next(row for row in read_csv(ROOT / "SOURCE_INDEX.csv") if row["file"] == mesh_name)
    actual = hashlib.sha256((folder / "mesh_original.mat").read_bytes()).hexdigest()
    if actual != record["sha256"]:
        raise ValueError("The published coarse mesh differs from its indexed checksum.")
    return record["original_sha256"]


def soil_samples(field):
    coordinates = np.vstack((field["soil_coordinates_m"], field["quadrature_coordinates_m"]))
    log_multipliers = np.r_[field["soil_log_multiplier"], field["quadrature_log_multiplier"]]
    return coordinates, log_multipliers


def compare_fields(folder, sample_id, mesh_digest):
    """Compare J7 radii and soil log multipliers at spatial locations shared by the meshes."""
    coarse = loadmat(folder / f"{sample_id}_original_field.mat", simplify_cells=True)
    fine = loadmat(folder / f"{sample_id}_refined_field.mat", simplify_cells=True)
    if coarse["mesh_sha256"] != mesh_digest:
        raise ValueError(f"The coarse field for {sample_id} refers to a different mesh.")
    if not np.array_equal(coarse["joint_radius_profile_m"], fine["joint_radius_profile_m"]):
        raise ValueError(f"The J7 radius profile changed on refinement for {sample_id}.")

    coarse_xy, coarse_values = soil_samples(coarse)
    fine_xy, fine_values = soil_samples(fine)
    coarse_by_xy = {tuple(xy): value for xy, value in zip(coarse_xy, coarse_values)}
    shared = [i for i, xy in enumerate(fine_xy) if tuple(xy) in coarse_by_xy]
    if len(fine_values) != 49665 or len(shared) != 45655:
        raise ValueError(f"Unexpected numbers of refined or shared soil locations for {sample_id}.")
    expected = [coarse_by_xy[tuple(fine_xy[i])] for i in shared]
    if not np.array_equal(fine_values[shared], expected):
        raise ValueError(f"Soil log multipliers changed at shared locations for {sample_id}.")


def rolling_one_day_peak(end_days, rates, reported_start):
    """Check the peak at every candidate boundary and at the reported start by direct integration."""
    edges = np.r_[0, end_days]
    widths = np.diff(edges)
    if not np.all(widths > 0) or end_days[-1] != HORIZON:
        raise ValueError("Leakage history must span five years in ascending time order.")
    # An interval edge can enter or leave the 24-hour window. Between those
    # events, the integral is linear in the starting time.
    starts = np.unique(np.r_[0.0, edges, edges - 1.0, HORIZON - 1])
    starts = starts[(starts >= 0) & (starts <= HORIZON - 1)]
    cumulative_volume = np.r_[0.0, np.cumsum(rates * widths)]
    window_volumes = (np.interp(starts + 1, edges, cumulative_volume)
                      - np.interp(starts, edges, cumulative_volume))
    peak = float(np.max(window_volumes))  # A one-day volume divided by one day.
    overlap = np.maximum(0, np.minimum(end_days, reported_start + 1) - np.maximum(edges[:-1], reported_start))
    if not abs(overlap @ rates - peak) < 1e-10:
        raise ValueError("The reported 24-hour window does not reproduce the maximum.")
    return peak


def verify(folder):
    folder = Path(folder)
    reported = read_csv(folder / "peaks.csv")
    threshold = float(reported[0]["threshold_L_m_d"])
    selected = nearest_cases(read_csv(folder / "selection.csv"), threshold)
    expected_pairs = {(sid, mesh) for sid in selected for mesh in ("original", "refined")}
    if len(reported) != len(expected_pairs) or {(r["case_id"], r["mesh"]) for r in reported} != expected_pairs:
        raise ValueError("Expected an original and a refined mesh result for each selected case.")
    mesh_digest = original_mesh_digest(folder)
    for sid in selected:
        compare_fields(folder, sid, mesh_digest)

    with np.load(folder / "rolling_rate_histories.npz", allow_pickle=False) as histories:
        for row in reported:
            case, mesh = row["case_id"], row["mesh"]
            if float(row["threshold_L_m_d"]) != threshold:
                raise ValueError(f"The threshold differs for {case} on the {mesh} mesh.")
            times, rates = histories[f"{case}_{mesh}"].T
            peak = rolling_one_day_peak(times, rates, float(row["window_start_day"]))
            if not abs(peak - float(row["peak_L_m_d"])) < 1e-10:
                raise ValueError(f"The integrated peak differs from the reported peak for {case} ({mesh}).")
            margin = 100 * (peak / threshold - 1)
            above = int(peak > threshold)
            if above != int(row["exceeds"]) or not abs(margin - float(row["margin_pct"])) < 1e-7:
                raise ValueError(f"Threshold classification or margin differs for {case} ({mesh}).")
            print(case, mesh, f"{peak:.9f}", above)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, nargs="?", default=DATA,
                        help="Directory containing peaks.csv and the saved mesh comparison.")
    verify(parser.parse_args().folder)
