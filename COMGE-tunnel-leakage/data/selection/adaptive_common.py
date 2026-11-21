"""Adaptive-only paths, immutable sources, and once guards; no physical submit."""

import os

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "4"
from pathlib import Path
import sys, json, time, traceback

A = Path(__file__).resolve().parent
X = A.parent
S = X.parent
BASE = S.parent
P = X / "first72_postprocess_20261011"
sys.path.insert(0, str(P))
import common72 as q

c = q.c
import numpy as np

sys.path.insert(0, str(BASE / "daily_peak_research_20261009"))
from fit_daily_peak import predict

KINDS = ["exponential", "matern32", "squared_exponential"]


def parameter_key(s):
    return json.dumps(
        [
            round(float(s["soil_cov"]), 12),
            round(float(s["joint_cov"]), 12),
            *[round(float(v), 12) for v in s["soil_theta_m"]],
            round(float(s["joint_theta_m"]), 12),
            s["soil_kernel"],
            s["joint_kernel"],
        ],
        separators=(",", ":"),
    )


def metadata(s):
    return np.array(
        [
            s["soil_cov"],
            s["joint_cov"],
            *np.log(s["soil_theta_m"]),
            np.log(s["joint_theta_m"]),
            *[
                float(s[key] == kind)
                for key in ["soil_kernel", "joint_kernel"]
                for kind in KINDS
            ],
        ]
    )


def verify():
    q.verify_sources()
    proof = c.read(A / "source_contracts.json")
    if not proof["passed"]:
        raise ValueError("Adaptive contracts missing")
    for p, h in proof["files"].items():
        if c.sha(p) != h:
            raise ValueError("Adaptive source/input changed " + p)
    review = q.rebind_proof(P / "independent_training584_review.json")
    if not review["passed"]:
        raise ValueError("584 training review missing")
    return proof


def run_once(name, source, function):
    verify()
    p = A / (name + "_intent.json")
    if p.exists():
        raise RuntimeError(name + " already submitted; no repeat")
    c.create(
        p,
        dict(
            createdUTC=c.utc(),
            source_sha256=c.sha(source),
            contracts_sha256=c.sha(A / "source_contracts.json"),
            physical_calls=0,
            validation_used=False,
        ),
    )
    start = time.perf_counter()
    cpu = time.process_time()
    try:
        function(start, cpu)
    except BaseException:
        c.create(
            A / (name + "_failure.json"),
            dict(
                createdUTC=c.utc(),
                source_sha256=c.sha(source),
                traceback=traceback.format_exc(),
                wall_s=time.perf_counter() - start,
                process_cpu_s=time.process_time() - cpu,
                physical_calls=0,
                automatic_retry=False,
            ),
        )
        raise


def percent_ranks(values):
    values = np.asarray(values, float)
    return np.array([np.mean(values <= v) for v in values])


def select_pairs(rows):
    selected = []
    comparisons = []
    for family in range(12):
        group = [r for r in rows if r["family_index"] == family]
        pairs = []
        for i in range(4):
            for j in range(i + 1, 4):
                l, r = group[i], group[j]
                with np.load(A / "fields" / l["candidate_file"]) as f:
                    g1 = np.array(f["soil_gaussian"])
                with np.load(A / "fields" / r["candidate_file"]) as f:
                    g2 = np.array(f["soil_gaussian"])
                dx = (np.array(l["features"]) - r["features"]) / np.array(
                    l["training_feature_scale"]
                )
                contrast = float(np.sqrt(np.mean(((g1 - g1.mean()) - (g2 - g2.mean())) ** 2)))
                pairs.append(
                    dict(
                        candidate_indices=[l["candidate_index"], r["candidate_index"]],
                        feature_distance28=float(np.linalg.norm(dx)),
                        centered_full_soil_gaussian_RMS_difference=contrast,
                        mean_priority=(l["priority_score"] + r["priority_score"]) / 2,
                    )
                )
        near = percent_ranks(-np.array([r["feature_distance28"] for r in pairs]))
        space = percent_ranks([r["centered_full_soil_gaussian_RMS_difference"] for r in pairs])
        for n, r in enumerate(pairs):
            r["pair_score"] = 0.6 * r["mean_priority"] + 0.25 * near[n] + 0.15 * space[n]
        winner = max(pairs, key=lambda r: (r["pair_score"], -min(r["candidate_indices"])))
        selected.extend(winner["candidate_indices"])
        comparisons.append(
            dict(
                family_index=family,
                pairs=pairs,
                selected_candidate_indices=winner["candidate_indices"],
            )
        )
    return selected, comparisons
