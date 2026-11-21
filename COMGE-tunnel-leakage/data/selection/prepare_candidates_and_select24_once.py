"""48 full independent candidate inputs, select24 without physical/validation data."""

from adaptive_common import *
import joblib

sys.path.insert(0, str(X / "input_sources"))
import uq_workflow as w
from field_sampler import draw, smooth_covariance_certificate, covariance
from hybrid_kl import sample_exponential_raw


def work(start, cpu):
    diagnostic = c.read(A / "family_diagnostic.json")
    if not diagnostic["passed"]:
        raise ValueError("Family diagnostic missing")
    (A / "fields").mkdir(exist_ok=False)
    geometry = dict(np.load(X / "input_sources/geometry_factors.npz"))
    points = np.vstack([geometry["nodes"], geometry["quadrature"]])
    certs = {}
    rows = []
    training = dict(np.load(P / "training584.npz"))
    newmodel = joblib.load(P / "peak584_rbf.joblib")
    oldmodel = joblib.load(S / "peak512/peak512_rbf.joblib")
    scale = newmodel["input_scaler"].scale_
    oldseeds = set(c.read(A / "prior_input_seed_catalog.json")["unique_seeds"])
    allseeds = []
    for family in diagnostic["selected_families"]:
        s = family["parameters"]
        name = family["scenario_id"]
        soil = (
            dict(
                passed=True,
                kind="exponential",
                method="Original exact irregular OU on complete original coordinates",
                continuum_2percent_RMS_qualified=False,
            )
            if s["soil_kernel"] == "exponential"
            else smooth_covariance_certificate(points, s["soil_theta_m"], s["soil_kernel"])
        )
        C = covariance(
            geometry["joint_coordinates"],
            geometry["joint_coordinates"],
            [s["joint_theta_m"]],
            s["joint_kernel"],
        )
        eig = np.linalg.eigvalsh((C + C.T) / 2)
        if eig.min() < -1e-12 * eig.max():
            raise ValueError("Indefinite joint covariance")
        certs[name] = dict(
            passed=True,
            soil=soil,
            joint=dict(
                passed=True,
                minimum_eigenvalue=float(eig.min()),
                maximum_diagonal_defect=float(np.max(abs(np.diag(C) - 1))),
                kind=s["joint_kernel"],
                dimension=64,
                PSD_factor_for_SE=True,
            ),
            scope="Original finite discrete sampler; not response/continuum accuracy",
        )
        for k in range(4):
            n = family["family_index"] * 4 + k
            seeds = [
                int(
                    np.random.SeedSequence([20261011, 12802, n, stream]).generate_state(
                        1, dtype=np.uint32
                    )[0]
                )
                for stream in range(2)
            ]
            allseeds += seeds
            if any(v in oldseeds for v in seeds) or len(set(allseeds)) != len(allseeds):
                raise ValueError("New seed collision")
            c.atomic(
                A / "progress.json",
                dict(
                    updatedUTC=c.utc(),
                    state="full_candidate_input_preparation",
                    completed=len(rows),
                    current_candidate=n,
                    total=48,
                    physical_calls=0,
                ),
            )
            g = draw(
                points, s["soil_theta_m"], s["soil_kernel"], seeds[0], sample_exponential_raw
            )
            j = draw(
                geometry["joint_coordinates"],
                [s["joint_theta_m"]],
                s["joint_kernel"],
                seeds[1],
                sample_exponential_raw,
            )
            f = w.physical_features(g, j, s["soil_cov"], s["joint_cov"], geometry, s)
            file = A / "fields" / f"adaptive_candidate_{n:03d}.npz"
            with file.open("xb") as stream:
                np.savez_compressed(
                    stream,
                    soil_gaussian=g,
                    joint_gaussian=j,
                    physical_gp_features=f,
                    seeds=np.array(seeds, dtype=np.uint32),
                )
            p584 = float(predict(newmodel, f[None, :])[0])
            p512 = float(predict(oldmodel, f[None, :])[0])
            distance = float(np.linalg.norm((training["X28"] - f) / scale, axis=1).min())
            row = dict(
                candidate_index=n,
                family_index=family["family_index"],
                scenario_id=name,
                parameters=s,
                seeds=seeds,
                seed_namespace=[20261011, 12802, n],
                candidate_file=file.name,
                field_sha256=c.sha(file),
                features=f.tolist(),
                training_feature_scale=scale.tolist(),
                prediction584=p584,
                prediction512=p512,
                model_disagreement_pct=abs(p584 - p512) / p584 * 100,
                nearest_training_feature_distance=distance,
            )
            c.create(file.with_suffix(".json"), row)
            rows.append(row)
            print(json.dumps(dict(candidate=n + 1, total=48, scenario=name)), flush=True)
    c.create(A / "covariance_certificates.json", certs)
    coverage = percent_ranks([r["nearest_training_feature_distance"] for r in rows])
    disagreement = percent_ranks([r["model_disagreement_pct"] for r in rows])
    tail = percent_ranks([r["prediction584"] for r in rows])
    for n, r in enumerate(rows):
        r["priority_score"] = float(0.4 * coverage[n] + 0.35 * disagreement[n] + 0.25 * tail[n])
    chosen, pairs = select_pairs(rows)
    byid = {r["candidate_index"]: r for r in rows}
    requests = []
    for n, index in enumerate(chosen):
        r = byid[index]
        s = r["parameters"]
        descriptor = [
            0,
            0,
            s["soil_cov"],
            s["joint_cov"],
            *s["soil_theta_m"],
            s["joint_theta_m"],
            KINDS.index(s["soil_kernel"]),
            KINDS.index(s["joint_kernel"]),
            *r["seeds"],
        ]
        requests.append(
            dict(
                sample_id=41072 + n,
                request_key=f"A128_ADAPT{n:02d}",
                role="training",
                stage="adaptive24",
                scenario_id=r["scenario_id"],
                parameters=s,
                seeds=r["seeds"],
                seed_namespace=r["seed_namespace"],
                candidate_file=r["candidate_file"],
                candidate_sha256=r["field_sha256"],
                candidate_index=index,
                descriptor=descriptor,
                selection_priority_score=r["priority_score"],
            )
        )
    c.create(A / "candidate_scores.json", rows)
    c.create(A / "requests24.json", requests)
    c.create(
        A / "adaptive_selection_proof.json",
        dict(
            createdUTC=c.utc(),
            passed=True,
            source_sha256=c.sha(__file__),
            contracts_sha256=c.sha(A / "source_contracts.json"),
            family_diagnostic_sha256=c.sha(A / "family_diagnostic.json"),
            candidate_scores_sha256=c.sha(A / "candidate_scores.json"),
            requests24_sha256=c.sha(A / "requests24.json"),
            certificates_sha256=c.sha(A / "covariance_certificates.json"),
            candidate_fields=48,
            selected_fields=24,
            selected_candidate_indices=chosen,
            pair_comparisons=pairs,
            feature_sufficiency_proven=False,
            exact_feature_equivalence_claimed=False,
            validation_used=False,
            physical_calls=0,
            full_field_generation_wall_s_once=time.perf_counter() - start,
            process_cpu_s=time.process_time() - cpu,
        ),
    )
    c.atomic(
        A / "progress.json",
        dict(
            updatedUTC=c.utc(),
            state="selection_completed_pending_independent_review",
            completed=48,
            selected24=24,
            physical_calls=0,
        ),
    )
    print(
        json.dumps(dict(passed=True, candidates=48, selected=24, physical_calls=0)), flush=True
    )


if __name__ == "__main__":
    run_once("prepare_candidates_select24", __file__, work)
