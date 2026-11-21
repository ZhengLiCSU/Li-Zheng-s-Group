"""Weighted Sobol-Nystrom KL with a covariance-consistent residual.

Adaptation of the CG Nystrom-KLE/orthogonal-residual idea to irregular points.
The raw residual uses exact separable exponential OU covariance instead of
a regular-grid FFT spectral approximation. A conditional cross-covariance
correction is essential when the approximate KL subspace is not invariant.
No rank truncation is applied to the remaining fluctuations.
"""

import os

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_key, "4")
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from scipy.linalg import cholesky, cho_solve
from scipy.stats import qmc


def kernel(x, y, lengths):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ell = np.asarray(lengths, float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[1] != y.shape[1] or ell.shape != (x.shape[1],):
        raise ValueError("Coordinate and correlation length dimensions disagree")
    if (
        not np.all(np.isfinite(x))
        or not np.all(np.isfinite(y))
        or not np.all(np.isfinite(ell))
        or np.any(ell <= 0)
    ):
        raise ValueError("Finite coordinates and positive lengths required")
    r = np.zeros((len(x), len(y)))
    for k in range(x.shape[1]):
        r += abs(x[:, None, k] - y[None, :, k]) / ell[k]
    return np.exp(-r)


def _ou_axis(noise, axis, ell):
    """Exact nonuniform stationary OU row; chunking prevents exp overflow."""
    out = np.empty(len(axis))
    out[0] = noise[0]
    first = 1
    while first < len(axis):
        last = int(np.searchsorted(axis, axis[first - 1] + 8 * ell, side="right"))
        last = max(last, first + 1)
        t = (axis[first:last] - axis[first - 1]) / ell
        b = np.sqrt(-np.expm1(-2 * np.diff(axis[first - 1 : last]) / ell))
        out[first:last] = np.exp(-t) * (
            out[first - 1] + np.cumsum(np.exp(t) * b * noise[first:last])
        )
        first = last
    return out


def sample_exponential_raw(points, lengths, seed):
    """Exact Gaussian covariance at irregular 1D/2D points, bounded memory.

    A separable nonuniform OU tensor is streamed one row at a time. Cost is
    nx*ny innovations; memory is O(nx+number of requested points), not nx*ny.
    Coordinates are not rounded or interpolated. The caller must consider
    tensor size when budgeting generation time.
    """
    p = np.asarray(points, float)
    ell = np.asarray(lengths, float)
    if p.ndim != 2 or p.shape[1] not in (1, 2) or ell.shape != (p.shape[1],):
        raise ValueError("Only one- and two-dimensional exponential fields supported")
    if (
        len(p) == 0
        or not np.all(np.isfinite(p))
        or not np.all(np.isfinite(ell))
        or np.any(ell <= 0)
    ):
        raise ValueError("Invalid field domain")
    rng = np.random.default_rng(seed)
    x, ix = np.unique(p[:, 0], return_inverse=True)
    if p.shape[1] == 1:
        return _ou_axis(rng.standard_normal(len(x)), x, ell[0])[ix]
    y, iy = np.unique(p[:, 1], return_inverse=True)
    order = np.argsort(iy, kind="stable")
    counts = np.bincount(iy, minlength=len(y))
    ends = np.cumsum(counts)
    starts = np.r_[0, ends[:-1]]
    prev = None
    answer = np.empty(len(p))
    for j in range(len(y)):
        innovation = _ou_axis(rng.standard_normal(len(x)), x, ell[0])
        if prev is None:
            row = innovation
        else:
            a = np.exp(-(y[j] - y[j - 1]) / ell[1])
            b = np.sqrt(-np.expm1(-2 * (y[j] - y[j - 1]) / ell[1]))
            row = a * prev + b * innovation
        select = order[starts[j] : ends[j]]
        answer[select] = row[ix[select]]
        prev = row
    return answer


def sobol_landmarks(bounds, count, seed, inside=None):
    bounds = np.asarray(bounds, float)
    d = bounds.shape[0]
    if bounds.shape != (d, 2) or np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError("Invalid bounds")
    sampler = qmc.Sobol(d, scramble=True, seed=seed)
    accepted = []
    n = 0
    for _ in range(64):
        z = qmc.scale(sampler.random(256), bounds[:, 0], bounds[:, 1])
        if inside is not None:
            z = z[np.asarray(inside(z), bool)]
        accepted.append(z)
        n += len(z)
        if n >= count:
            return np.vstack(accepted)[:count]
    raise ValueError("Insufficient landmarks inside the actual soil domain")


@dataclass
class HybridKL:
    anchors: np.ndarray
    weights: np.ndarray
    lengths: np.ndarray
    U: np.ndarray
    eigenvalues: np.ndarray
    report: dict

    @classmethod
    def build(cls, anchors, weights, lengths, landmarks, head_modes=64):
        x = np.asarray(anchors, float)
        w = np.asarray(weights, float).reshape(-1)
        ell = np.asarray(lengths, float)
        if len(w) != len(x) or np.any(w <= 0) or not np.all(np.isfinite(w)):
            raise ValueError("Positive domain weights required")
        w = w / sum(w)
        sw = np.sqrt(w)
        z = np.asarray(landmarks, float)
        K = kernel(z, z, ell) / len(z)
        lam, V = eigh(K, driver="evd")
        order = np.argsort(lam)[::-1]
        lam, V = lam[order], V[:, order]
        if not isinstance(head_modes, int) or head_modes < 1 or len(z) < 1:
            raise ValueError("Positive integer head count and landmarks required")
        m = min(head_modes, len(z), len(x))
        if lam[m - 1] <= max(lam[0] * 1e-13, 0):
            raise ValueError("Landmark basis is numerically rank deficient")
        F = kernel(x, z, ell) @ V[:, :m] / np.sqrt(len(z) * lam[:m])
        U, R = np.linalg.qr(sw[:, None] * F, mode="reduced")
        if np.min(abs(np.diag(R))) < 1e-12:
            raise ValueError("Transferred Nystrom basis lost rank")
        BU = np.empty_like(U)
        for start in range(0, len(x), 256):
            BU[start : start + 256] = sw[start : start + 256, None] * (
                kernel(x[start : start + 256], x, ell) @ (sw[:, None] * U)
            )
        A = U.T @ BU
        asym = np.linalg.norm(A - A.T) / max(np.linalg.norm(A), np.finfo(float).tiny)
        if asym > 1e-10:
            raise ArithmeticError("Ritz covariance is not symmetric")
        r, Q = eigh((A + A.T) / 2)
        order = np.argsort(r)[::-1]
        r, Q = r[order], Q[:, order]
        if r[-1] <= 0:
            raise ArithmeticError("Nonpositive retained covariance")
        U = U @ Q
        BU = BU @ Q
        D = BU - U * r
        report = {
            "method": "Sobol Nystrom + weighted Ritz + conditional orthogonal residual",
            "head_modes": m,
            "landmarks": len(z),
            "anchors": len(x),
            "head_trace_fraction": float(sum(r)),
            "ritz_asymmetry": float(asym),
            "orthogonality_defect": float(np.linalg.norm(U.T @ U - np.eye(m))),
            "subspace_cross_covariance_norm": float(np.linalg.norm(D)),
            "head_only_gaussian_RMS": float(np.sqrt(max(0, 1 - sum(r)))),
            "residual_included": True,
            "raw_covariance": "Exact separable exponential at requested points",
            "continuum_FE_quadrature_accuracy_qualified": False,
        }
        return cls(x, w, ell, U, r, report)

    def evaluation_factor(self, points):
        """Conditional mean factor C_ea W^1/2 U / sqrt(lambda)."""
        x = np.asarray(points, float)
        factor = np.empty((len(x), len(self.eigenvalues)))
        right = np.sqrt(self.weights)[:, None] * self.U / np.sqrt(self.eigenvalues)
        for start in range(0, len(x), 256):
            factor[start : start + 256] = (
                kernel(x[start : start + 256], self.anchors, self.lengths) @ right
            )
        return factor

    def sample(self, points, head_seed, residual_seed, head_coefficients=None, factor=None):
        p = np.asarray(points, float)
        xi = (
            np.random.default_rng(head_seed).standard_normal(len(self.eigenvalues))
            if head_coefficients is None
            else np.asarray(head_coefficients, float).reshape(-1)
        )
        if len(xi) != len(self.eigenvalues) or not np.all(np.isfinite(xi)):
            raise ValueError("Invalid head coefficients")
        raw = sample_exponential_raw(np.vstack([self.anchors, p]), self.lengths, residual_seed)
        anchor_raw = raw[: len(self.anchors)]
        value_raw = raw[len(self.anchors) :]
        old_head = (self.U.T @ (np.sqrt(self.weights) * anchor_raw)) / np.sqrt(self.eigenvalues)
        F = self.evaluation_factor(p) if factor is None else np.asarray(factor, float)
        if F.shape != (len(p), len(xi)):
            raise ValueError("Evaluation factor shape mismatch")
        return value_raw + F @ (xi - old_head), xi

    def anchor_sample(self, head_seed, residual_seed, chol=None, factor=None):
        """Cheap candidate features; the same anchors are retained in a full run."""
        L = (
            cholesky(kernel(self.anchors, self.anchors, self.lengths), lower=True)
            if chol is None
            else chol
        )
        raw = L @ np.random.default_rng(residual_seed).standard_normal(len(self.anchors))
        xi = np.random.default_rng(head_seed).standard_normal(len(self.eigenvalues))
        old = (self.U.T @ (np.sqrt(self.weights) * raw)) / np.sqrt(self.eigenvalues)
        F = self.evaluation_factor(self.anchors) if factor is None else factor
        return raw + F @ (xi - old), xi, raw

    def extend_anchor_sample(
        self, points, xi, raw_anchor, extension_seed, chol=None, factor=None
    ):
        """Matheron conditioning extends an existing candidate without reseeding it.

        A candidate's anchor residual is sampled cheaply with dense Cholesky.
        An independent exact OU draw is then conditioned on those SAME raw
        anchor values. This prevents the common error of changing a residual
        realization when more quadrature coordinates are requested.
        """
        p = np.asarray(points, float)
        xi = np.asarray(xi, float)
        a = np.asarray(raw_anchor, float)
        if xi.shape != self.eigenvalues.shape or a.shape != (len(self.anchors),):
            raise ValueError("Anchor state shape mismatch")
        L = (
            cholesky(kernel(self.anchors, self.anchors, self.lengths), lower=True)
            if chol is None
            else chol
        )
        draw = sample_exponential_raw(
            np.vstack([self.anchors, p]), self.lengths, extension_seed
        )
        beta = cho_solve((L, True), a - draw[: len(a)])
        values = draw[len(a) :].copy()
        for start in range(0, len(p), 256):
            values[start : start + 256] += (
                kernel(p[start : start + 256], self.anchors, self.lengths) @ beta
            )
        old = (self.U.T @ (np.sqrt(self.weights) * a)) / np.sqrt(self.eigenvalues)
        F = self.evaluation_factor(p) if factor is None else np.asarray(factor, float)
        if F.shape != (len(p), len(xi)):
            raise ValueError("Factor shape mismatch")
        return values + F @ (xi - old)

    def covariance_certificate(self, points, physical_cov=0.5):
        """Deterministic algebra defect of the finite-point Gaussian covariance.

        This verifies covariance, NOT pathwise error to a different sampler,
        continuum integration error, or Monte Carlo estimator precision.
        """
        p = np.asarray(points, float)
        C = kernel(p, p, self.lengths)
        F = self.evaluation_factor(p)
        A = np.empty((len(self.U), len(self.eigenvalues)))
        right = np.sqrt(self.weights)[:, None] * self.U
        for start in range(0, len(self.U), 256):
            A[start : start + 256] = np.sqrt(self.weights[start : start + 256, None]) * (
                kernel(self.anchors[start : start + 256], self.anchors, self.lengths) @ right
            )
        H = (self.U.T @ A) / np.sqrt(self.eigenvalues[:, None] * self.eigenvalues[None, :])
        approximate = C + F @ (H - np.eye(len(H))) @ F.T
        s2 = np.log1p(physical_cov**2)
        ref = np.expm1(s2 * C)
        appr = np.expm1(s2 * approximate)
        return {
            "gaussian_covariance_relative_Frobenius_defect": float(
                np.linalg.norm(approximate - C) / np.linalg.norm(C)
            ),
            "lognormal_covariance_relative_Frobenius_defect": float(
                np.linalg.norm(appr - ref) / np.linalg.norm(ref)
            ),
            "maximum_point_variance_defect": float(max(abs(np.diag(approximate) - 1))),
            "scope": "Prescribed finite-point covariance; not a full continuum RMS/FE quadrature certificate",
            "RMS_2percent_total_qualified": False,
        }

    def save(self, file):
        import json

        np.savez_compressed(
            file,
            anchors=self.anchors,
            weights=self.weights,
            lengths=self.lengths,
            U=self.U,
            eigenvalues=self.eigenvalues,
            report=json.dumps(self.report),
        )

    @classmethod
    def load(cls, file):
        import json

        with np.load(file) as d:
            return cls(
                d["anchors"],
                d["weights"],
                d["lengths"],
                d["U"],
                d["eigenvalues"],
                json.loads(str(d["report"])),
            )


def lognormal_field(gaussian, arithmetic_mean, cov):
    if (
        not np.isfinite(arithmetic_mean)
        or arithmetic_mean <= 0
        or not np.isfinite(cov)
        or cov <= 0
    ):
        raise ValueError("Positive finite mean and COV required")
    a = np.log1p(cov * cov)
    log_value = np.log(arithmetic_mean) + np.sqrt(a) * np.asarray(gaussian) - a / 2
    values = np.exp(log_value)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ArithmeticError("Field outside floating point range")
    return values, log_value
