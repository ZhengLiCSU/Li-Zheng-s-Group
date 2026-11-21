"""Frozen multi-family input materialization. Features in labels are request
descriptors, NOT surrogate regressors. Physical GP regressors are separately
derived from each immutable full input field after staging.
"""

import os

for k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(k, "4")
from pathlib import Path
import hashlib, json, time
import numpy as np
from scipy.io import savemat
from field_sampler import REVISION, draw, physical_fields, smooth_covariance_certificate
from hybrid_kl import sample_exponential_raw


def digest(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def write(p, x):
    with Path(p).open("x", encoding="utf-8") as f:
        json.dump(x, f, ensure_ascii=False, indent=2, allow_nan=False)


def load_cache(cache):
    cache = Path(cache)
    m = json.loads((cache / "cache_manifest.json").read_text(encoding="utf-8"))
    for p, h in m["files"].items():
        if digest(cache / p) != h:
            raise ValueError("Cache changed " + p)
    for p, h in m["sources_sha256"].items():
        if digest(p) != h:
            raise ValueError("Generator source changed " + p)
    return m


def physical_features(full, j, soil_cov, joint_cov, d, scenario):
    """Return the 28 Kriging inputs for one soil and segment-joint realization.

    The soil field is ordered by node, then by quadrature point. The 64 J7
    radii follow the flow path through the lining. Both fields are supplied
    as standard Gaussian values before the lognormal transformation.
    """
    from scipy.special import logsumexp

    soil_log_multiplier, joint_radius_m, _ = physical_fields(full, j, soil_cov, joint_cov)
    log_radius = np.log(joint_radius_m / 1e-6)

    # The inverse second and fourth powers weight small radii along J7.
    joint_inputs = [
        float(-(logsumexp(-2 * log_radius) - np.log(64))),
        float(-(logsumexp(-4 * log_radius) - np.log(64))),
        *log_radius.reshape(4, 16).mean(axis=1).tolist(),
        float(log_radius.std()),
    ]

    node_ids = d["residual_witness_ids"]
    coordinates = d["nodes"][node_ids]
    soil_values = soil_log_multiplier[node_ids]
    # Coordinates on the upper edge belong to the last of the three bins.
    relative_position = (coordinates - coordinates.min(axis=0)) / np.ptp(coordinates, axis=0)
    region = np.minimum((3 * relative_position).astype(int), 2)
    regional_means = []
    for x_bin in range(3):
        for y_bin in range(3):
            in_region = (region[:, 0] == x_bin) & (region[:, 1] == y_bin)
            regional_means.append(float(soil_values[in_region].mean()))

    # Input order also fixes the columns of the frozen training matrix.
    inputs = [
        *joint_inputs,
        *regional_means,
        float(soil_values.std()),
        soil_cov,
        joint_cov,
        *np.log(scenario["soil_theta_m"]).tolist(),
        float(np.log(scenario["joint_theta_m"])),
    ]
    kernel_order = ("exponential", "matern32", "squared_exponential")
    for property_name in ("soil_kernel", "joint_kernel"):
        inputs.extend(float(scenario[property_name] == kernel) for kernel in kernel_order)

    if len(inputs) != 28 or not np.all(np.isfinite(inputs)):
        raise ValueError("A field realization must supply 28 finite Kriging inputs.")
    return np.asarray(inputs)


def materialize(cache, sample_id, destination):
    start = time.perf_counter()
    cache = Path(cache)
    m = load_cache(cache)
    r = m["requests"][str(sample_id)]
    s = m["scenarios"][r["scenario_id"]]
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    d = dict(np.load(cache / "geometry_factors.npz"))
    points = np.vstack([d["nodes"], d["quadrature"]])
    joint = d["joint_coordinates"]
    g = draw(points, s["soil_theta_m"], s["soil_kernel"], r["seeds"][0], sample_exponential_raw)
    j = draw(
        joint, [s["joint_theta_m"]], s["joint_kernel"], r["seeds"][1], sample_exponential_raw
    )
    soil_log, radius, rmean = physical_fields(g, j, s["soil_cov"], s["joint_cov"])
    profile = np.full((8, 64), 1e-6)
    profile[6] = radius
    phy = physical_features(g, j, s["soil_cov"], s["joint_cov"], d, s)
    f = np.asarray(r["descriptor"], float)
    payload = dict(
        soil_coordinates_m=d["nodes"],
        soil_log_multiplier=soil_log[: len(d["nodes"]), None],
        quadrature_coordinates_m=d["quadrature"],
        quadrature_log_multiplier=soil_log[len(d["nodes"]) :, None],
        joint_radius_profile_m=profile,
        features=f[None, :],
        sample_id=sample_id,
        generator_revision=REVISION,
        cache_manifest_sha256=digest(cache / "cache_manifest.json"),
        mesh_sha256=m["mesh_sha256"],
        seed_vector=np.asarray(r["seeds"], dtype=np.uint32),
        joint_radial_coordinates_m=joint,
        soil_Ksat_mean_m_s=1e-7,
        physical_gp_features=phy[None, :],
        scenario_id=r["scenario_id"],
    )
    if not np.all(np.isfinite(soil_log)) or np.max(abs(soil_log)) > 100:
        raise ArithmeticError("Unrepresentable input field")
    savemat(destination / "field_input.mat", payload, do_compression=True)
    np.savez_compressed(
        destination / "realization.npz",
        features=f,
        physical_gp_features=phy,
        soil_gaussian=g,
        joint_gaussian=j,
        seeds=r["seeds"],
    )
    report = dict(
        revision=REVISION,
        sample_id=sample_id,
        seeds=r["seeds"],
        cache=str(cache.resolve()),
        cache_manifest_sha256=digest(cache / "cache_manifest.json"),
        scenario_id=r["scenario_id"],
        scenario=s,
        request_role=r["role"],
        descriptor_is_not_GP_input=True,
        physical_feature_count=28,
        field_input_sha256=digest(destination / "field_input.mat"),
        realization_sha256=digest(destination / "realization.npz"),
        materialization_wall_s=time.perf_counter() - start,
        soil_node_count=len(d["nodes"]),
        soil_SCV_count=len(d["quadrature"]),
        expected_radius_arithmetic_mean_m=rmean,
        fixed_pointwise_radius_fourth_moment_m4=1e-24 * (1 + 0.25**2) ** 6,
        effective_harmonic_conductance_mean_fixed=False,
        features_include_residual=True,
        full_continuum_2percent_RMS_qualified=False,
        covariance_certificate=m["covariance_certificates"][r["scenario_id"]],
        solver_invoked=False,
    )
    write(destination / "realization_manifest.json", report)
    return report


def read_labels(files, cache_hash):
    rows = []
    seen = set()
    for f in files:
        p = Path(f)
        d = json.loads(p.read_text(encoding="utf-8"))
        if not all(
            d.get(k) is True
            for k in ("completeFiveYearHistory", "physicalAuditPassed", "sourceHashVerified")
        ):
            raise ValueError("Incomplete physical label")
        if (
            d["cache_manifest_sha256"] != cache_hash
            or digest(p.parent / "engineering_five_year_result.mat") != d["result_sha256"]
        ):
            raise ValueError("Result/cache binding changed")
        i = int(d["sample_id"])
        y = float(d["Q5_L_per_m"])
        if i in seen or not np.isfinite(y) or y < 0:
            raise ValueError("Duplicate or invalid response")
        seen.add(i)
        rows.append((i, y))
    return rows
