"""Reconstruct the field features, Kriging predictions, and leakage statistics in the paper.

Run from a local copy of the study with Python, NumPy and SciPy:
    python -B verify_review.py --output-dir review_check

The output directory contains the held-out predictions, the density histogram
for Fig. 12(d), and verification.json with the numerical comparisons.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

# Apply defaults before importing NumPy, while respecting the caller's settings.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
REVISION = DATA / "revision"
THRESHOLD = 0.17827224170071532  # L m⁻¹ d⁻¹, full precision used in the parameter study
REPLAY_TOLERANCE = 1e-10        # L m⁻¹ d⁻¹
REPLICATIONS = 8
HORIZON = 1825                  # days
sys.path.insert(0, str(ROOT / "code/original"))
from uq_workflow import physical_features


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, header, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def check_manifest():
    hashes = read_json(ROOT / "manifest.json")
    for name, expected in hashes.items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Checksum mismatch for {name}")
    return len(hashes)


def reconstruct_features():
    """Recompute the 28 model inputs for the stored full-field example."""
    example = DATA / "example"
    parameters = read_json(example / "adaptive_candidate_000.json")["parameters"]
    with np.load(example / "adaptive_candidate_000.npz", allow_pickle=False) as fields, \
         np.load(example / "geometry.npz", allow_pickle=False) as geometry:
        features = physical_features(
            fields["soil_gaussian"], fields["joint_gaussian"],
            parameters["soil_cov"], parameters["joint_cov"], dict(geometry), parameters,
        )
        max_error = float(np.max(np.abs(features - fields["physical_gp_features"])))
    if not max_error <= 1e-12:
        raise ValueError(f"Stored and reconstructed field features differ by {max_error:g}.")
    return features, max_error


def predict(inputs, model):
    """Return five-year peak leakage rates (L m⁻¹ d⁻¹) from the frozen Kriging fit."""
    scaled_inputs = (inputs - model["input_mean"]) / model["input_scale"]
    lengths = model["length_scale"]
    distances_squared = cdist(scaled_inputs / lengths, model["X_train_scaled"] / lengths, "sqeuclidean")
    kernel = float(model["signal_variance"]) * np.exp(-0.5 * distances_squared)
    # The training target was standardized log leakage; reverse both transforms.
    log_leakage = float(model["target_mean"]) + float(model["target_scale"]) * (kernel @ model["alpha"])
    return np.exp(log_leakage)


def check_training_folds(training_ids):
    """Every training case must appear in a test fold once under each split scheme."""
    schemes = read_json(DATA / "model/training608_evaluation.json")["splits"]
    for scheme, folds in schemes.items():
        held_out = []
        for number, fold in enumerate(folds, start=1):
            train, test = set(fold["train_ids"]), set(fold["test_ids"])
            if train & test or (train | test) != training_ids:
                raise ValueError(f"Overlapping or missing cases in {scheme}, fold {number}")
            if scheme == "parameter_group8" and set(fold["train_groups"]) & set(fold["test_groups"]):
                raise ValueError(f"A parameter group occurs on both sides of fold {number}.")
            held_out.extend(fold["test_ids"])
        if sorted(held_out) != sorted(training_ids):
            raise ValueError(f"The test folds in {scheme} do not partition the training cases.")


def replay_held_out(model, output):
    with np.load(DATA / "validation/validation80_dataset.npz", allow_pickle=False) as test, \
         np.load(DATA / "model/training608.npz", allow_pickle=False) as training:
        training_ids = set(training["ids"])
        if len(training_ids) != 608 or len(set(test["ids"])) != 80 or len(test["ids"]) != 80:
            raise ValueError("Expected 608 distinct training cases and 80 distinct held-out cases.")
        if training_ids.intersection(test["ids"]):
            raise ValueError("A held-out case also occurs in the training set.")
        check_training_folds(training_ids)
        simulated = test["y"]
        predicted = predict(test["X28"], model)
        replay_error = float(np.max(np.abs(predicted - test["prediction608"])))
        if not replay_error <= REPLAY_TOLERANCE:
            raise ValueError(f"Maximum replay error {replay_error:g} exceeds {REPLAY_TOLERANCE:g}.")
        write_csv(output / "validation_replayed.csv",
                  ["sample_id", "cohort", "physical_L_m-1_d-1", "replayed_L_m-1_d-1"],
                  zip(test["ids"], test["groups"], simulated, predicted))

    subsets = {"all80": np.ones(len(simulated), dtype=bool),
               "near_threshold_within20percent": np.abs(simulated / THRESHOLD - 1) <= 0.2}
    metrics = {}
    for name, mask in subsets.items():
        true, estimate = simulated[mask], predicted[mask]
        metrics[name] = {
            "n": int(mask.sum()),
            "physical_exceedances": int(np.sum(true > THRESHOLD)),
            "relative_L2_percent": float(100 * np.linalg.norm(estimate - true) / np.linalg.norm(true)),
            "maximum_absolute_relative_error_percent": float(100 * np.max(np.abs(estimate / true - 1))),
            "classification_disagreements": int(np.sum((estimate > THRESHOLD) != (true > THRESHOLD))),
        }
    return metrics, replay_error


def compare_sample_sizes():
    """Maximum probability (percentage points), mean (%) and 95th percentile (%) changes."""
    rows = read_csv(REVISION / "mc_stability_cases.csv")
    estimates = {(row["case_id"], int(row["N"])): row for row in rows}
    cases = sorted({row["case_id"] for row in rows})
    sizes = (8192, 16384, 32768, 65536, 131072)
    if len(cases) != 4181 or len(rows) != len(cases) * len(sizes) or len(rows) != len(estimates):
        raise ValueError("Expected five sample sizes for each of the 4,181 parameter settings.")
    changes = []
    for case in cases:
        for size in sizes:
            row = estimates[case, size]
            count = int(row["exceedance_count"])
            if not 0 <= count <= size or not abs(float(row["probability"]) - count / size) < 1e-15:
                raise ValueError(f"Exceedance count and probability disagree for {case}, N={size}.")
        lower, final = estimates[case, 65536], estimates[case, 131072]
        changes.append([100 * abs(float(lower["probability"]) - float(final["probability"])),
                        100 * abs(float(lower["mean"]) / float(final["mean"]) - 1),
                        100 * abs(float(lower["p95"]) / float(final["p95"]) - 1)])
    return np.max(changes, axis=0).tolist(), len(cases)


def check_paired_probabilities():
    """Recompute contrasts and their standard errors across paired Gaussian draws."""
    pairs = read_csv(REVISION / "mc_stability_paired.csv")
    with np.load(REVISION / "mc_stability_replicates.npz", allow_pickle=False) as saved:
        by_case = {str(case): i for i, case in enumerate(saved["case_ids"])}
        by_size = {int(n) * REPLICATIONS: i for i, n in enumerate(saved["n_per_rep"])}
        probabilities = saved["prob_mean_p95"][..., 0]  # case × sample size × replication
        for row in pairs:
            n = by_size[int(row["N"])]
            left, right = row["left_case"], row["right_case"]
            delta = probabilities[by_case[right], n] - probabilities[by_case[left], n]
            mean, se = delta.mean(), delta.std(ddof=1) / np.sqrt(REPLICATIONS)
            if not abs(mean - float(row["delta_p"])) < 1e-15 or not abs(se - float(row["delta_p_SE"])) < 1e-14:
                raise ValueError(f"Paired probability contrast mismatch for {left} and {right}, N={row['N']}.")
    return len(pairs)


def integrate_mesh_history(case):
    """Integrate piecewise-constant interval rates; times are days and volumes are L m⁻¹."""
    name = case["case"]
    with np.load(REVISION / f"{name}_response.npz", allow_pickle=False) as history:
        edges, rate = history["interval_edges"], history["rate"]
        if len(edges) != len(rate) + 1 or edges[0] != 0 or not abs(edges[-1] - HORIZON) < 1e-8:
            raise ValueError(f"Invalid five-year time interval in {name}.")
        widths = np.diff(edges)
        if not np.all(widths > 0):
            raise ValueError(f"Unordered time steps in {name}.")
        volume = np.r_[0, np.cumsum(widths * rate)]
        # A rolling one-day integral can change slope only when either end reaches an interval edge.
        starts = np.clip(np.r_[edges, edges - 1], 0, HORIZON - 1)
        peak = float(np.max(np.interp(starts + 1, edges, volume) - np.interp(starts, edges, volume)))
        daily = np.diff(np.interp(np.arange(HORIZON + 1), edges, volume))
        year_five = float(np.sum(daily[-365:]))
        np.testing.assert_allclose(daily, history["daily"], atol=1e-10, rtol=0, err_msg=name)
        np.testing.assert_allclose([peak, year_five], [case["peak_L_m_d"], case["annual_volumes_L_m"][-1]],
                                   atol=1e-9, rtol=0, err_msg=name)
    return peak, year_five, daily


def check_mesh_comparisons():
    result = read_json(REVISION / "matched_physics_results.json")
    histories = {case["case"]: integrate_mesh_history(case) for case in result["cases"]}
    for pair in result["contrasts"]:
        field = pair["field"]
        peak0, vol0, daily0 = histories[field + "_M0"]
        peak1, vol1, daily1 = histories[field + "_M1"]
        np.testing.assert_allclose([peak0, peak1, vol0, vol1],
                                   [pair["M0_peak"], pair["M1_peak"], pair["M0_year5_volume"], pair["M1_year5_volume"]],
                                   atol=1e-9, rtol=0, err_msg=field)
        daily_error = 100 * np.sum(np.abs(daily0 - daily1)) / np.sum(daily1)
        np.testing.assert_allclose(daily_error, pair["five_year_daily_L1_percent"], atol=1e-8, rtol=0, err_msg=field)
        if not abs(100 * (peak0 / peak1 - 1) - pair["peak_difference_percent"]) < 1e-10:
            raise ValueError(f"Five-year peak difference disagrees for {field}.")
        if not abs(100 * (vol0 / vol1 - 1) - pair["volume_difference_percent"]) < 1e-10:
            raise ValueError(f"Year-five volume difference disagrees for {field}.")
    return len(histories)


def reproduce_figure12(output):
    """Recreate the 55-bin density histogram and exceedance probability in Fig. 12(d)."""
    with np.load(DATA / "figure12/PC01739.npz", allow_pickle=False) as source:
        peaks = source["predictions"].reshape(-1)
    case = read_json(DATA / "figure12/PC01739_completed.json")
    if len(peaks) != 131072 or not abs(peaks.mean() - case["mean"]) < 1e-14:
        raise ValueError("Figure 12(d) sample count or mean differs from the recorded result.")
    if not abs(float(np.quantile(peaks, 0.95)) - case["p95"]) < 1e-14:
        raise ValueError("Figure 12(d) 95th percentile differs from the recorded result.")
    density, edges = np.histogram(peaks, bins=55, density=True)
    write_csv(output / "figure12d_histogram.csv",
              ["left_L_m-1_d-1", "right_L_m-1_d-1", "density_m_d_L-1"],
              zip(edges[:-1], edges[1:], density))
    return case["parameters"], float(np.mean(peaks > THRESHOLD))


def verify(output):
    count = check_manifest()
    example, feature_error = reconstruct_features()
    with np.load(DATA / "model/frozen_predictor_arrays.npz", allow_pickle=False) as fitted:
        model = dict(fitted)
    validation, replay_error = replay_held_out(model, output)
    max_change, settings = compare_sample_sizes()
    pairs = check_paired_probabilities()
    histories = check_mesh_comparisons()
    fig_parameters, fig_probability = reproduce_figure12(output)
    return {
        "passed": True, "files_checked": count, "feature_reconstruction_max_abs": feature_error,
        "prediction_replay_max_abs_L_m-1_d-1": replay_error, "predictor_tolerance_abs": REPLAY_TOLERANCE,
        "example_prediction_L_m-1_d-1": float(predict(example[None, :], model)[0]), "validation": validation,
        "MC_max_change_probability_pp_mean_percent_p95_percent": max_change, "MC_settings": settings,
        "paired_contrasts_recomputed": pairs, "CV_splits_checked": True,
        "physical_histories_reintegrated": histories, "figure12d_parameters": fig_parameters,
        "figure12d_probability": fig_probability, "new_solver_calls": 0, "new_model_fits": 0, "new_random_fields": 0,
        "scope": "Reconstruction of saved inputs, frozen predictions and statistics, not a new physical or probability calibration experiment.",
    }


def main(out):
    target = Path(out).resolve()
    if target.exists():
        raise FileExistsError(f"Choose a fresh output directory: {target} already exists.")
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=target.parent) as scratch:
        output = Path(scratch) / "results"
        output.mkdir()
        summary = json.dumps(verify(output), indent=2)
        (output / "verification.json").write_text(summary, encoding="utf-8")
        output.rename(target)
    print(summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("review_check"),
                        help="New directory for verification.json and the two CSV outputs.")
    main(parser.parse_args().output_dir)
