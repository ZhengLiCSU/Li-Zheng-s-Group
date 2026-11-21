"""Train-only nested CV for daily peaks, then freeze predictions on existing32 inputs."""

import os

for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[k] = "4"
from pathlib import Path
import csv, hashlib, json, time, warnings, sys
from datetime import datetime, timezone
import numpy as np
import joblib
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
from prepare_daily_peak import HERE, BASE, VALID, read, write, sha, utc

sys.path.insert(0, str(BASE / "surrogate_accuracy_20261009"))
from feature_map import make_parameters, transform_features

warnings.filterwarnings("ignore", category=ConvergenceWarning)


def metrics(y, p):
    rel = abs(p - y) / y * 100
    top = y >= np.quantile(y, 0.9)
    return dict(
        relative_L2_pct=float(np.linalg.norm(p - y) / np.linalg.norm(y) * 100),
        RMSE_L_m_d=float(np.sqrt(np.mean((p - y) ** 2))),
        MAE_L_m_d=float(np.mean(abs(p - y))),
        max_relative_pct=float(rel.max()),
        p95_relative_pct=float(np.quantile(rel, 0.95)),
        bias_L_m_d=float(np.mean(p - y)),
        top10pct_relative_L2_pct=float(
            np.linalg.norm((p - y)[top]) / np.linalg.norm(y[top]) * 100
        ),
        top10pct_underprediction_max_pct=float(
            np.max(np.maximum(y[top] - p[top], 0) / y[top]) * 100
        ),
    )


def optimizer(obj, theta, bounds):
    opt = minimize(
        obj,
        theta,
        jac=True,
        bounds=bounds,
        method="L-BFGS-B",
        options={"maxiter": 100, "ftol": 1e-8},
    )
    return opt.x, opt.fun


def fit(X, y, kernel):
    scaler = StandardScaler().fit(X)
    z = scaler.transform(X)
    t = np.log(y)
    mu = float(t.mean())
    sd = float(t.std())
    correlation = (
        RBF(np.ones(z.shape[1]) * 3, (0.05, 1e4))
        if kernel == "rbf"
        else Matern(np.ones(z.shape[1]) * 3, (0.05, 1e4), nu=2.5)
    )
    gp = GaussianProcessRegressor(
        kernel=ConstantKernel(1.0, (1e-3, 1e3)) * correlation + WhiteKernel(0.01, (1e-6, 0.5)),
        normalize_y=False,
        optimizer=optimizer,
        alpha=1e-9,
        random_state=20261009,
    ).fit(z, (t - mu) / sd)
    gp.optimizer = None
    return dict(input_scaler=scaler, gp=gp, target_mean=mu, target_scale=sd, kernel=kernel)


def predict(packet, X):
    t = packet["gp"].predict(packet["input_scaler"].transform(X))
    return np.exp(packet["target_mean"] + packet["target_scale"] * t)


def main():
    start = time.perf_counter()
    data = dict(np.load(HERE / "training_daily_peaks.npz"))
    X = data["physical4"]
    fold = data["fold"]
    ids = data["sample_ids"]
    provenance = read(HERE / "training_data_provenance.json")
    assert (
        provenance["passed"]
        and sha(HERE / "training_daily_peaks.npz") == provenance["dataset_sha256"]
    )
    plan = read(HERE / "research_plan.json")
    assert plan["training_count"] == 384
    assert not (HERE / "frozen_daily_peak_models.joblib").exists()
    report = {}
    packets = {}
    oofs = []
    modelrows = []
    final_cfg = []
    for ti, target in enumerate(data["targets"]):
        target = str(target)
        y = data["y"][:, ti]
        p = np.full(384, np.nan)
        mobility = np.full(384, np.nan)
        ratio = np.full(384, np.nan)
        records = []
        for f in range(1, 9):
            tr = np.where(fold != f)[0]
            te = np.where(fold == f)[0]
            inner = KFold(3, shuffle=True, random_state=20261009 + f)
            scores = []
            inner_records = []
            for kernel in ("rbf", "matern"):
                ip = np.full(len(tr), np.nan)
                for it, iv in inner.split(tr):
                    ip[iv] = predict(fit(X[tr[it]], y[tr[it]], kernel), X[tr[iv]])
                scores.append(float(np.mean((ip - y[tr]) ** 2)))
                inner_records.append(dict(kernel=kernel, MSE=scores[-1]))
            selected = ("rbf", "matern")[int(np.argmin(scores))]
            packet = fit(X[tr], y[tr], selected)
            p[te] = predict(packet, X[te])
            lr = LinearRegression().fit(X[tr, :1], np.log(y[tr]))
            mobility[te] = np.exp(lr.predict(X[te, :1]))
            ratio[te] = data["Q5"][te] * np.mean(y[tr] / data["Q5"][tr])
            records.append(
                dict(
                    fold=f,
                    train_ids=ids[tr].tolist(),
                    test_ids=ids[te].tolist(),
                    inner_splits=[
                        dict(train_ids=ids[tr[it]].tolist(), test_ids=ids[tr[iv]].tolist())
                        for it, iv in inner.split(tr)
                    ],
                    candidates=inner_records,
                    selected=selected,
                    kernel=str(packet["gp"].kernel_),
                    metrics=metrics(y[te], p[te]),
                )
            )
            print(
                f"{target} fold{f}/8 {selected} relL2={metrics(y[te], p[te])['relative_L2_pct']:.4f}%, elapsed={time.perf_counter() - start:.1f}s",
                flush=True,
            )
        assert np.all(np.isfinite(p))
        # Final selection uses only384 training CV, never the existing32 targets.
        scores = []
        for kernel in ("rbf", "matern"):
            ip = np.full(384, np.nan)
            for it, iv in KFold(3, shuffle=True, random_state=20261009).split(X):
                ip[iv] = predict(fit(X[it], y[it], kernel), X[iv])
            scores.append(float(np.mean((ip - y) ** 2)))
        selected = ("rbf", "matern")[int(np.argmin(scores))]
        packet = fit(X, y, selected)
        packet.update(
            target=target,
            training_ids=ids.copy(),
            physical_changes=False,
            feature_name="physical4",
        )
        packets[target] = packet
        oofs.append(p)
        final_cfg.append(selected)
        report[target] = dict(
            nested8fold=metrics(y, p),
            input_only_mobility_powerlaw=metrics(y, mobility),
            annualQ5_to_peak_ratio_diagnostic=metrics(y, ratio),
            ratio_diagnostic_uses_true_Q5_not_input_only=True,
            peak_to_Q5_ratio=dict(
                min=float(np.min(y / data["Q5"])),
                max=float(np.max(y / data["Q5"])),
                mean=float(np.mean(y / data["Q5"])),
                CV_pct=float(np.std(y / data["Q5"]) / np.mean(y / data["Q5"]) * 100),
            ),
            response_range=dict(
                min=float(y.min()),
                max=float(y.max()),
                mean=float(y.mean()),
                CV_pct=float(y.std() / y.mean() * 100),
            ),
            final_kernel=selected,
            final3fold_MSE=scores,
            final_kernel_parameters=str(packet["gp"].kernel_),
            fold_records=records,
        )
        modelrows.extend(
            dict(
                sample_id=int(ids[i]),
                fold=int(fold[i]),
                target=target,
                true_L_m_d=float(y[i]),
                OOF_L_m_d=float(p[i]),
                mobility_powerlaw_OOF_L_m_d=float(mobility[i]),
                Q5_ratio_diagnostic_OOF_L_m_d=float(ratio[i]),
            )
            for i in range(384)
        )
    manifest = read(VALID / "campaign_manifest.json")
    pool = np.load(VALID / "design/candidate_pool.npz")
    vids = np.asarray(manifest["validation_ids"])
    assert len(vids) == 32 and not set(vids) & set(ids)
    cache = Path(manifest["frozen_cache"])
    params = make_parameters(np.load(cache / "geometry_factors.npz"))
    VX = transform_features(pool["features"][vids], params, 4)
    assert np.allclose(X, transform_features(data["X"], params, 4), rtol=1e-13, atol=1e-13)
    vp = np.column_stack([predict(packets[str(t)], VX) for t in data["targets"]])
    np.savez_compressed(
        HERE / "frozen_existing32_daily_predictions.npz",
        sample_ids=vids,
        targets=data["targets"],
        predicted=vp,
        physical4=VX,
    )
    joblib.dump(packets, HERE / "frozen_daily_peak_models.joblib")
    np.savez_compressed(
        HERE / "nested8fold_daily_predictions.npz",
        sample_ids=ids,
        targets=data["targets"],
        fold=fold,
        true=data["y"],
        predicted=np.column_stack(oofs),
    )
    with (HERE / "nested8fold_daily_predictions.csv").open(
        "x", encoding="utf-8-sig", newline=""
    ) as f:
        w = csv.DictWriter(f, list(modelrows[0]))
        w.writeheader()
        w.writerows(modelrows)
    write(
        HERE / "training_evaluation.json",
        dict(
            createdUTC=utc(),
            passed=True,
            source_sha256=sha(__file__),
            training_dataset_sha256=sha(HERE / "training_daily_peaks.npz"),
            targets=report,
            models_sha256=sha(HERE / "frozen_daily_peak_models.joblib"),
            predictions32_sha256=sha(HERE / "frozen_existing32_daily_predictions.npz"),
            existing32_daily_responses_read=False,
            existing32_used_for_tuning=False,
            all_transforms_fit_inside_training_partitions=True,
            known_limit="Adaptive training design; randomized8fold is interpolation assessment, not deployment/tail/reliability certification. Existing32 previously inspected for Q5; new peak check remains exploratory.",
            wall_s=time.perf_counter() - start,
        ),
    )
    print("DAILY_MODELS_AND_32_PREDICTIONS_FROZEN", flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=4):
        main()
