"""Prespecified log RBF ARD28 with within-fold scaling and parameter-family CV."""

import os

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "4"
from common24_post import *
import numpy as np, joblib
from sklearn.model_selection import KFold, GroupKFold

sys.path.insert(0, str(BASE / "daily_peak_research_20261009"))
from fit_daily_peak import fit, predict, metrics


def work(start, cpu):
    proof = rebind_proof(P / "dataset608_proof.json")
    if c.sha(P / "training608.npz") != proof["dataset_sha256"]:
        raise ValueError("608 dataset changed")
    d = dict(np.load(P / "training608.npz"))
    X28 = d["X28"]
    y = d["y"]
    ids = d["ids"]
    groups = d["parameter_groups"]
    preds = {}
    reports = {}
    splits = {}
    for name, splitter in [
        ("random8", KFold(8, shuffle=True, random_state=20261011)),
        ("parameter_group8", GroupKFold(8)),
    ]:
        out = np.full(len(y), np.nan)
        folds = []
        for k, (tr, te) in enumerate(splitter.split(X28, y, groups)):
            packet = fit(X28[tr], y[tr], "rbf")
            out[te] = predict(packet, X28[te])
            folds.append(
                dict(
                    fold=k,
                    train_ids=ids[tr].tolist(),
                    test_ids=ids[te].tolist(),
                    train_groups=sorted(set(groups[tr])),
                    test_groups=sorted(set(groups[te])),
                    input_scaler_mean=packet["input_scaler"].mean_.tolist(),
                    input_scaler_scale=packet["input_scaler"].scale_.tolist(),
                    target_mean=packet["target_mean"],
                    target_scale=packet["target_scale"],
                    kernel=str(packet["gp"].kernel_),
                )
            )
            c.atomic(
                P / "fit_progress.json",
                dict(
                    updatedUTC=c.utc(),
                    scheme=name,
                    completed_fold=k + 1,
                    total_folds=8,
                    physical_calls=0,
                ),
            )
            print(
                json.dumps(dict(scheme=name, fold=k + 1, wall_s=time.perf_counter() - start)),
                flush=True,
            )
        if not np.isfinite(out).all():
            raise ValueError("Missing OOF predictions")
        preds[name] = out
        splits[name] = folds
        reports[name] = {
            key: metrics(y[mask], out[mask])
            for key, mask in [
                ("all608", np.ones(len(y), bool)),
                ("new96", ids >= 41000),
                ("first72", (ids >= 41000) & (ids < 41072)),
                ("adaptive24", ids >= 41072),
                ("old512", ids < 41000),
            ]
        }
    final = fit(X28, y, "rbf")
    joblib.dump(final, P / "peak608_rbf.joblib")
    with (P / "cv608_predictions.npz").open("xb") as f:
        np.savez_compressed(f, ids=ids, y=y, **preds)
    c.create(
        P / "training608_evaluation.json",
        dict(
            createdUTC=c.utc(),
            passed=True,
            source_sha256=c.sha(__file__),
            dataset_proof_sha256=c.sha(P / "dataset608_proof.json"),
            source_contracts_sha256=c.sha(P / "source_contracts.json"),
            model_sha256=c.sha(P / "peak608_rbf.joblib"),
            cv_sha256=c.sha(P / "cv608_predictions.npz"),
            count=608,
            metrics=reports,
            splits=splits,
            full_kernel=str(final["gp"].kernel_),
            wall_s=time.perf_counter() - start,
            process_cpu_s=time.process_time() - cpu,
            validation_used=False,
            physical_calls=0,
            feature_sufficiency_proven=False,
            GP_variance_calibrated=False,
            independent_validation_accuracy_claimed=False,
        ),
    )
    print(
        json.dumps(
            dict(passed=True, count=608, metrics=reports, wall_s=time.perf_counter() - start)
        ),
        flush=True,
    )


if __name__ == "__main__":
    run_once("fit_peak608", __file__, work)
