"""Gaussian marginal of frozen discrete feature inputs; no full-field/solver claim."""

import os

for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "4"
import sys, json, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, eigh
from scipy.fft import fft, ifft, next_fast_len
from scipy.special import logsumexp

S = Path(__file__).resolve().parent
B = S.parent / "uq_batches/b96_20261009_rf_reliability_stage01_v3"
M = json.loads((B / "campaign_manifest.json").read_text(encoding="utf-8"))
sys.path.insert(0, M["frozen_uq"])
import field_sampler as fs
import uq_workflow as uq


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, x):
    with Path(p).open("x", encoding="utf-8") as f:
        json.dump(x, f, ensure_ascii=False, indent=2, allow_nan=False)


def geometry():
    return dict(np.load(Path(M["frozen_cache"]) / "geometry_factors.npz"))


def scenarios():
    return json.loads((B / "design/approved_study_plan.json").read_text(encoding="utf-8"))[
        "all_scenarios"
    ]


def axis_cov(indices, n, step, ell, kind):
    """Actual grid covariance, including the implemented finite FFT clipping."""
    delta = np.abs(indices[:, None] - indices[None, :])
    if kind == "squared_exponential":
        N = next_fast_len(max(2 * (int(n) - 1), int(np.ceil(16 * ell / step))))
        h = np.minimum(np.arange(N), N - np.arange(N)) * step
        ev = fft(fs.rho(h, ell, kind)).real
        if ev.min() < -1e-9:
            raise ArithmeticError("Original embedding rejection")
        actual = ifft(np.maximum(ev, 0)).real
        return actual[delta], dict(embedding_size=N, min_unclipped_eigenvalue=float(ev.min()))
    # Actual stationary state recursion, retaining finite arithmetic Q/Cholesky.
    t = step / ell
    A = np.exp(-t) * np.array([[1 + t, t], [-t, 1 - t]])
    Q = np.eye(2) - A @ A.T
    L = cholesky((Q + Q.T) / 2, lower=True)
    V = [np.eye(2)]
    for k in range(int(indices.max())):
        V.append(A @ V[-1] @ A.T + L @ L.T)
    powers = {int(k): np.linalg.matrix_power(A, int(k)) for k in np.unique(delta)}
    C = np.empty(delta.shape)
    for a, i in enumerate(indices):
        for b, j in enumerate(indices):
            C[a, b] = (powers[abs(int(i) - int(j))] @ V[min(int(i), int(j))])[0, 0]
    return C, dict(
        max_stationary_state_roundoff=float(max(np.max(abs(v - np.eye(2))) for v in V))
    )


def soil_covariance(d, s):
    points = np.vstack([d["nodes"], d["quadrature"]])
    p = d["nodes"][d["residual_witness_ids"]]
    kind = s["soil_kernel"]
    theta = s["soil_theta_m"]
    if kind == "exponential":
        return fs.covariance(p, p, theta, kind), dict(kind=kind, full_geometry_grid=False)
    origin, n, step, ell = fs.grid_spec(points, theta, kind)
    i, a = fs.interpolation(p, origin, n, step)
    offsets = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    corners = (i[:, None, :] + offsets).reshape(-1, 2)
    weights = np.c_[
        (1 - a[:, 0]) * (1 - a[:, 1]),
        (1 - a[:, 0]) * a[:, 1],
        a[:, 0] * (1 - a[:, 1]),
        a[:, 0] * a[:, 1],
    ]
    axes = []
    proof = []
    for k in range(2):
        ids, inv = np.unique(corners[:, k], return_inverse=True)
        C, q = axis_cov(ids, n[k], step[k], ell[k], kind)
        axes.append(C[inv[:, None], inv[None, :]])
        proof.append(q)
    C = (axes[0] * axes[1]).reshape(len(p), 4, len(p), 4)
    C = np.einsum("ia,iajb,jb->ij", weights, C, weights)
    variance = np.prod(
        (1 - a) ** 2
        + a * a
        + 2 * a * (1 - a) * np.array([fs.rho(step[k], ell[k], kind) for k in range(2)]),
        axis=1,
    )
    C /= np.sqrt(variance[:, None] * variance[None, :])
    return C, dict(
        kind=kind,
        full_geometry_grid=True,
        origin=origin.tolist(),
        n=n.tolist(),
        step=step.tolist(),
        axes=proof,
        extra_variance_renormalization=False,
    )


def factor(C):
    values, vectors = eigh((C + C.T) / 2)
    if values.min() < -1e-12 * values.max():
        raise ArithmeticError("Indefinite marginal covariance")
    F = vectors * np.sqrt(np.maximum(values, 0))[None, :]
    defect = float(np.linalg.norm(F @ F.T - C) / np.linalg.norm(C))
    if defect > 1e-12:
        raise ArithmeticError("Marginal reconstruction defect")
    return F, defect


def features(g, j, d, s):
    """Batch equivalent of the original28 feature function, on its dependencies."""
    logk, rad, _ = fs.physical_fields(g, j, s["soil_cov"], s["joint_cov"])
    lr = np.log(rad / 1e-6)
    xy = d["nodes"][d["residual_witness_ids"]]
    zone = np.minimum(((xy - xy.min(0)) / np.ptp(xy, axis=0) * 3).astype(int), 2)
    cols = [
        -(logsumexp(-2 * lr, axis=1) - np.log(64)),
        -(logsumexp(-4 * lr, axis=1) - np.log(64)),
    ]
    cols.extend(lr.reshape(-1, 4, 16).mean(2).T)
    cols.append(lr.std(1))
    cols.extend(
        logk[:, (zone[:, 0] == a) & (zone[:, 1] == b)].mean(1)
        for a in range(3)
        for b in range(3)
    )
    cols.append(logk.std(1))
    fixed = [
        s["soil_cov"],
        s["joint_cov"],
        *np.log(s["soil_theta_m"]),
        np.log(s["joint_theta_m"]),
    ]
    fixed.extend(
        float(s[key] == kind)
        for key in ("soil_kernel", "joint_kernel")
        for kind in ("exponential", "matern32", "squared_exponential")
    )
    cols.extend(np.full(len(g), v) for v in fixed)
    X = np.asarray(cols).T
    if X.shape != (len(g), 28) or not np.all(np.isfinite(X)):
        raise ArithmeticError("Invalid28 features")
    return X
