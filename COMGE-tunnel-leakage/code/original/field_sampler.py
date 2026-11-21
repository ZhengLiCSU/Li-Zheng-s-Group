"""New statistical-family sampler; no solver, response or old registry access.

Exponential fields use the archived exact irregular-coordinate OU sampler.
Smooth separable fields use stationary grid draws and explicit bilinear
evaluation. The finite-point interpolation covariance is certified per family;
this is not a continuum/response/reliability convergence certificate.
"""

import os

for k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(k, "4")
import numpy as np
from scipy.fft import rfft2, irfft2, fft, next_fast_len
from scipy.linalg import cholesky, eigh

REVISION = "RF_reliability_equal_SOF_fixed_local_r4_20261009_v1"


def lengths(kind, theta):
    theta = np.asarray(theta, float)
    if np.any(theta <= 0):
        raise ValueError("Positive Gaussian-field SOF required")
    return (
        theta
        / {"exponential": 2.0, "matern32": 4.0, "squared_exponential": np.sqrt(np.pi)}[kind]
    )


def rho(h, ell, kind):
    t = np.abs(np.asarray(h)) / ell
    if kind == "exponential":
        return np.exp(-t)
    if kind == "matern32":
        return (1 + t) * np.exp(-t)
    if kind == "squared_exponential":
        return np.exp(-t * t)
    raise ValueError(kind)


def covariance(x, y, theta, kind):
    ell = lengths(kind, theta)
    out = np.ones((len(x), len(y)))
    for k in range(len(ell)):
        out *= rho(x[:, None, k] - y[None, :, k], ell[k], kind)
    return out


def grid_spec(points, theta, kind):
    ell = lengths(kind, theta)
    origin = points.min(0)
    span = np.ptp(points, axis=0)
    n = np.maximum(2, np.ceil(span / (ell / 64)).astype(int) + 1)
    step = np.where(span > 0, span / (n - 1), ell / 64)
    return origin, n, step, ell


def interpolation(points, origin, n, step):
    z = (points - origin) / step
    lo = np.minimum(np.floor(z).astype(int), n - 2)
    lo = np.maximum(lo, 0)
    a = z - lo
    if points.shape[1] == 1:
        return lo[:, 0], a[:, 0]
    return lo, a


def smooth_covariance_certificate(points, theta, kind):
    origin, n, step, ell = grid_spec(points, theta, kind)
    select = np.unique(np.linspace(0, len(points) - 1, min(96, len(points))).astype(int))
    p = points[select]
    if p.shape[1] == 1:
        i, a = interpolation(p, origin, n, step)
        corners = origin[0] + (i[:, None] + np.arange(2)) * step[0]
        weights = np.c_[1 - a, a]
        x = corners.reshape(-1, 1)
    else:
        i, a = interpolation(p, origin, n, step)
        off = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
        x = (origin + (i[:, None, :] + off) * step).reshape(-1, 2)
        weights = np.c_[
            (1 - a[:, 0]) * (1 - a[:, 1]),
            (1 - a[:, 0]) * a[:, 1],
            a[:, 0] * (1 - a[:, 1]),
            a[:, 0] * a[:, 1],
        ]
    width = weights.shape[1]
    C = covariance(x, x, theta, kind).reshape(len(p), width, len(p), width)
    approx = np.einsum("ia,iajb,jb->ij", weights, C, weights)
    target = covariance(p, p, theta, kind)
    variance = float(np.max(abs(np.diag(approx) - 1)))
    approx = approx / np.sqrt(np.diag(approx)[:, None] * np.diag(approx)[None, :])
    defect = float(np.linalg.norm(approx - target) / np.linalg.norm(target))
    # Smooth kernels have a bounded mixed finite-difference interpolation
    # error; the explicit witness covariance is supplemented by step bounds.
    assert np.all(step / ell <= 1 / 64 + 1e-12) and variance < 0.001 and defect < 0.002
    return dict(
        passed=True,
        kind=kind,
        grid_shape=n.tolist(),
        grid_step_over_ell=(step / ell).tolist(),
        witness_count=len(p),
        max_witness_variance_defect=variance,
        witness_relative_covariance_defect=defect,
        pointwise_variance_normalized_to_one=True,
        finite_point_interpolation_used=True,
        continuum_2percent_RMS_qualified=False,
        response_convergence_qualified=False,
    )


def matern_grid(n, step, ell, rng):
    shape = tuple(n)
    state = np.empty((2, *shape))
    # Stationary normalized [value, derivative/lambda] state, covariance I.
    # First x-axis: columns are independent stationary Matérn3/2 paths.
    t = step[0] / ell[0]
    A = np.exp(-t) * np.array([[1 + t, t], [-t, 1 - t]])
    Q = np.eye(2) - A @ A.T
    L = cholesky((Q + Q.T) / 2, lower=True)
    prev = rng.standard_normal((2, 2, n[1]))
    state[:, 0, :] = prev[:, 0, :]
    for i in range(1, n[0]):
        prev = np.einsum("ab,pby->pay", A, prev) + np.einsum(
            "ab,pby->pay", L, rng.standard_normal(prev.shape)
        )
        state[:, i, :] = prev[:, 0, :]
    # Two independent spatial x paths drive the y-state innovations.
    t = step[1] / ell[1]
    A = np.exp(-t) * np.array([[1 + t, t], [-t, 1 - t]])
    Q = np.eye(2) - A @ A.T
    L = cholesky((Q + Q.T) / 2, lower=True)
    for j in range(1, n[1]):
        state[:, :, j] = A @ state[:, :, j - 1] + L @ state[:, :, j]
    return state[0]


def gaussian_grid(n, step, ell, rng):
    size = [
        next_fast_len(max(2 * (int(n[k]) - 1), int(np.ceil(16 * ell[k] / step[k]))))
        for k in range(2)
    ]
    eigen = []
    for k, N in enumerate(size):
        h = np.minimum(np.arange(N), N - np.arange(N)) * step[k]
        v = fft(rho(h, ell[k], "squared_exponential")).real
        if v.min() < -1e-9:
            raise ArithmeticError(
                "Circulant embedding not positive; no silent covariance repair"
            )
        eigen.append(np.maximum(v, 0))
    Z = rng.standard_normal(size)
    F = rfft2(Z, workers=4)
    del Z
    F *= np.sqrt(eigen[0])[:, None] * np.sqrt(eigen[1][: size[1] // 2 + 1])[None, :]
    out = irfft2(F, s=size, workers=4)
    return out[: n[0], : n[1]].copy()


def draw(points, theta, kind, seed, exact_exponential=None):
    p = np.asarray(points, float)
    ell = lengths(kind, theta)
    rng = np.random.default_rng(seed)
    if p.ndim != 2 or not np.all(np.isfinite(p)):
        raise ValueError("Finite coordinates required")
    if p.shape[1] == 1:
        C = covariance(p, p, theta, kind)
        if kind != "squared_exponential":
            return cholesky(C, lower=True) @ rng.standard_normal(len(p))
        values, vectors = eigh((C + C.T) / 2)
        if values.min() < -1e-12 * values.max():
            raise ArithmeticError("Joint covariance genuinely indefinite")
        factor = vectors * np.sqrt(np.maximum(values, 0))[None, :]
        if np.linalg.norm(factor @ factor.T - C) / np.linalg.norm(C) > 1e-12:
            raise ArithmeticError("Joint covariance reconstruction failed")
        return factor @ rng.standard_normal(len(p))
    if kind == "exponential":
        if exact_exponential is None:
            raise ValueError("Archived exact OU sampler required")
        return exact_exponential(p, ell, seed)
    origin, n, step, ell = grid_spec(p, theta, kind)
    G = (
        matern_grid(n, step, ell, rng)
        if kind == "matern32"
        else gaussian_grid(n, step, ell, rng)
    )
    i, a = interpolation(p, origin, n, step)
    x, y = i.T
    ax, ay = a.T
    values = (
        (1 - ax) * (1 - ay) * G[x, y]
        + (1 - ax) * ay * G[x, y + 1]
        + ax * (1 - ay) * G[x + 1, y]
        + ax * ay * G[x + 1, y + 1]
    )
    variance = ((1 - ax) ** 2 + ax**2 + 2 * ax * (1 - ax) * rho(step[0], ell[0], kind)) * (
        (1 - ay) ** 2 + ay**2 + 2 * ay * (1 - ay) * rho(step[1], ell[1], kind)
    )
    return values / np.sqrt(variance)


def physical_fields(g, j, soil_cov, joint_cov):
    s2 = np.log1p(soil_cov**2)
    soil_log = np.sqrt(s2) * g - s2 / 2
    # Keep ensemble pointwise E[r^4] equal to the original cv=.25 family.
    # This does not fix each realization's harmonic effective conductance.
    j2 = np.log1p(joint_cov**2)
    rmean = 1e-6 * ((1 + 0.25**2) / (1 + joint_cov**2)) ** 1.5
    radius = np.exp(np.log(rmean) + np.sqrt(j2) * j - j2 / 2)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0):
        raise ArithmeticError("Invalid positive radius")
    return soil_log, radius, rmean
