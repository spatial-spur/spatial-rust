"""Small-sample reference calculations using explicit matrices and integration."""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.spatial.distance import cdist


def distances(coordinates, geographic=False):
    if not geographic:
        return cdist(coordinates, coordinates)
    lon, lat = np.asarray(coordinates).T
    lon = np.remainder(lon + 180, 360) - 180
    at_pole = np.abs(lat) == 90
    lon[at_pole] = 0
    dlon = np.remainder(np.subtract.outer(lon, lon) + 180, 360) - 180
    dlat = np.subtract.outer(lat, lat)
    cos_lat = np.cos(np.deg2rad(lat))
    cos_lat[at_pole] = 0
    hav = np.sin(np.deg2rad(dlat / 2)) ** 2
    hav += np.outer(cos_lat, cos_lat) * np.sin(np.deg2rad(dlon / 2)) ** 2
    return np.arcsin(np.sqrt(np.clip(hav, 0, 1))) / np.pi


def center(n):
    return np.eye(n) - 1 / n


def brownian_gls(y, d):
    """Whiten the complete positive Brownian covariance subspace."""
    n = d.shape[0]
    centering = center(n)
    values, vectors = np.linalg.eigh(-0.5 * centering @ (d / d.max()) @ centering)
    keep = values > 1e-10
    transform = (vectors[:, keep] / np.sqrt(values[keep])) @ vectors[:, keep].T
    return transform @ (centering @ y)


def distance_kernels(d, avc):
    n = len(d)
    pairs = d[np.tril_indices(n, -1)]
    floor = np.mean(pairs == 0)
    if not floor < avc < 1:
        raise ValueError("Unattainable average correlation")

    def decay(rho):
        def f(rate):
            return np.exp(-rate * pairs).mean() - rho

        end = 1.0
        while f(end) > 0:
            end *= 2
        return brentq(f, 0, end, xtol=1e-13)

    c0 = decay(avc)
    cmax = decay(floor + min(1e-5 * (1 - floor), 0.1 * (avc - floor)))
    steps = max(2, int(np.ceil(np.log(cmax / c0) / np.log(1.05))))
    scales = np.unique(np.r_[c0 * 1.05 ** np.arange(steps), cmax])
    return [np.exp(-rate * d) for rate in scales] + [np.asarray(d == 0, dtype=float)]


def rejection(omega, statistic):
    """Gaussian ratio tail via the characteristic function of a quadratic form."""
    values, vectors = np.linalg.eigh(omega)
    root = vectors * np.sqrt(np.maximum(values, 0))
    q = len(values) - 1
    restriction = np.diag(np.r_[1.0, np.full(q, -(statistic**2) / q)])
    eigen = np.linalg.eigvalsh(root.T @ restriction @ root)
    eigen /= np.max(np.abs(eigen))
    eigen = eigen[np.abs(eigen) > 1e-14]
    if np.all(eigen <= 0):
        return 0.0
    if np.all(eigen >= 0):
        return 1.0
    if len(eigen) == 2:
        return 2 * np.arctan(np.sqrt(-eigen[1] / eigen[0])) / np.pi

    def integrand(t):
        x = 2 * t * eigen
        phase = 0.5 * np.arctan(x).sum()
        amplitude = np.exp(-0.25 * np.log1p(x * x).sum())
        return amplitude * np.sin(phase) / t

    value, _ = quad(integrand, 0, np.inf, epsabs=2e-10, epsrel=2e-10)
    return 0.5 + value / np.pi


def critical(omegas, alpha=0.05):
    def objective(t):
        return max(rejection(sigma, t) for sigma in omegas) - alpha

    upper = 2.0
    while objective(upper) > 0:
        upper *= 2
    return brentq(objective, 1e-7, upper, xtol=3e-8)


def regression(y, x, z=None):
    xhat = x if z is None else z @ np.linalg.lstsq(z, x, rcond=None)[0]
    b = np.linalg.lstsq(xhat, y, rcond=None)[0]
    residual = y - x @ b
    direction = xhat @ np.linalg.inv(xhat.T @ xhat)
    residual_transpose = np.eye(len(x)) - direction @ x.T
    return b, residual, direction, residual_transpose


def conditional_weights(w, direction, projection, membership):
    scales = np.linalg.norm(membership * direction[:, None], axis=0)
    counts = membership.sum(axis=0)
    signs = np.divide(
        direction,
        membership @ scales,
        out=np.zeros_like(direction),
        where=(membership @ scales) != 0,
    )
    signs[(membership @ scales == 0) & (membership @ counts == 1)] = 1
    row_weights = membership @ w
    numerator = row_weights[:, :1] * direction[:, None]
    denominator = projection @ (row_weights[:, 1:] * direction[:, None])
    return membership.T @ (signs[:, None] * np.column_stack((numerator, denominator)))


def scpc_at_dimension(y, x, d, q, avc, *, z=None, membership=None):
    """Fixed-q dense oracle for clustered and unclustered OLS or 2SLS."""
    if membership is None:
        membership = np.eye(len(y))
    kernels = distance_kernels(d, avc)
    p = center(len(d))
    eigenvectors = np.linalg.eigh(p @ kernels[0] @ p)[1]
    w = np.column_stack((np.full(len(d), 1 / np.sqrt(len(d))), eigenvectors[:, -q:]))
    estimate, residual, direction, transpose = regression(y, x, z)
    scores = membership.T @ (direction * residual[:, None])
    se = np.linalg.norm(w[:, 1:].T @ scores, axis=0) * np.sqrt(len(d) / q)
    location = [w.T @ kernel @ w for kernel in kernels]
    location_cutoff = critical(location)
    t = estimate / se
    pvalues, cutoffs, endpoints = [], [], []
    for j, stat in enumerate(t):
        wc = conditional_weights(w, direction[:, j], transpose, membership)
        conditional = [wc.T @ kernel @ wc for kernel in kernels]
        families = location + conditional
        cutoffs.append(critical(families))
        endpoints.append(critical([location[-1], conditional[-1]]))
        pvalues.append(max(rejection(omega, stat) for omega in families))
    return dict(
        coef=estimate,
        se=se,
        t=t,
        p=pvalues,
        critical=cutoffs,
        endpoint=endpoints,
        location_cutoff=location_cutoff,
        weights=w,
    )


def angular_log_density(covariance, z):
    """Scale-invariant Gaussian angular density, omitting constants in the dimension."""
    q = covariance.shape[0]
    _, logdet = np.linalg.slogdet(covariance)
    squared = np.sum(z * np.linalg.solve(covariance, z), axis=0)
    return -0.5 * (logdet + q * np.log(squared))


def rank_pvalue(observed, simulated):
    return (1 + np.count_nonzero(simulated >= observed)) / (len(simulated) + 1)


def i1_calibration(y, d, q, rate, draws):
    """Dense levels I(1) reference given the alternative rate selected by tuning."""
    d = d / d.max()
    p = center(len(d))
    brownian = -0.5 * p @ d @ p
    values, vectors = np.linalg.eigh(brownian)
    basis = vectors[:, -q:][:, ::-1]
    for col in basis.T:
        absolute = np.abs(col)
        anchor = col[np.flatnonzero(absolute >= (1 - 1e-12) * absolute.max())[0]]
        col *= np.sign(anchor)
    alternative = basis.T @ np.exp(-rate * d) @ basis
    null = np.diag(values[-q:][::-1])
    l0 = np.linalg.cholesky(null)
    la = np.linalg.cholesky(alternative)

    def ratio(z):
        return np.sum(np.linalg.solve(l0, z) ** 2, axis=0) / np.sum(
            np.linalg.solve(la, z) ** 2, axis=0
        )

    statistic = ratio(basis.T @ (y - y.mean()))
    simulated = ratio(l0 @ draws)
    ranked = np.sort(simulated)[::-1]
    ranks = [(len(simulated) + 1) // quantile for quantile in [100, 20, 10]]
    return (
        statistic,
        rank_pvalue(statistic, simulated),
        [ranked[rank - 1] if rank else np.inf for rank in ranks],
    )
