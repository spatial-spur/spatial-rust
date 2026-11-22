"""Compare finite calibration grids with independent Gaussian calculations."""

import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.spatial.distance import pdist, squareform

from scpc.utils.spatial import (
    get_c0_from_avc,
    get_cv,
    get_oms,
    get_w_with_ties,
    kernel_scale_grid,
    set_oms_wfin,
)
from spur.utils.inference import get_ha_param_i0
from spur.utils.matrix import get_r


def scpc_grid_case():
    sites = np.array(
        [
            0.09304880031893936,
            2.7447794507671857,
            2.752585040503658,
            3.007933601595235,
            3.0812773934737323,
            3.97084917611615,
        ]
    )
    x = np.array(
        [
            -0.018250213817659094,
            -1.581390328927324,
            -0.3472816224011446,
            1.4374529717030404,
            0.36442888490765174,
            0.14504030853543554,
        ]
    )
    d = squareform(pdist(sites[:, None]))
    distances = pdist(sites[:, None])
    c0 = get_c0_from_avc(distances, 0.5)
    cmax = get_c0_from_avc(distances, 1e-5)
    w, _ = get_w_with_ties(d, c0, 1)
    qx, _ = np.linalg.qr(np.column_stack((np.ones(len(x)), x)))
    residual = x * w[:, 1]
    residual -= qx @ (qx.T @ residual)
    wx = np.column_stack((np.abs(x) * w[:, 0], np.sign(x) * residual))

    def cauchy(weights, rate):
        om = weights.T @ np.exp(-rate * d) @ weights
        return om[0, 1] / om[1, 1], np.sqrt(np.linalg.det(om)) / om[1, 1]

    def cutoff(weights, rate):
        loc, scale = cauchy(weights, rate)
        v = scale / math.tan(math.pi * 0.05)
        return v + math.sqrt(v * v + loc * loc + scale * scale)

    rates = kernel_scale_grid(c0, cmax, 1.2)
    cv = max(cutoff(weights, rate) for weights in (w, wx) for rate in rates)
    off_grid = 1.1505417694486828
    loc, scale = cauchy(wx, off_grid)
    size = 1 - (math.atan((cv - loc) / scale) + math.atan((cv + loc) / scale)) / math.pi
    assert c0 < off_grid < cmax and size > 0.05
    print(f"SCPC fixed q=1 at decay {off_grid:.6f}: cutoff {cv:.6f}, size {size:.8%}")

    setup = set_oms_wfin(sites[:, None], 0.5, False, method="exact")
    w = setup.wfin
    q = w.shape[1] - 1
    residual = x[:, None] * w[:, 1:]
    residual -= qx @ (qx.T @ residual)
    wx = np.column_stack((np.abs(x) * w[:, 0], np.sign(x)[:, None] * residual))
    cv = max(setup.cvfin, get_cv(get_oms(d, c0, cmax, wx, 1.2), q, 0.05))
    off_grid = 0.8033512400507625
    omega = wx.T @ np.exp(-off_grid * d) @ wx
    chol = np.linalg.cholesky(0.5 * (omega + omega.T))
    eig = np.linalg.eigvalsh(chol.T @ np.diag([1] + [-(cv**2) / q] * q) @ chol)
    ratios = -eig[:-1] / eig[-1]
    size = (
        2
        / np.pi
        * quad(
            lambda theta: np.prod(1 + ratios / math.cos(theta) ** 2) ** -0.5,
            0,
            np.pi / 2,
        )[0]
    )
    assert rates[0] < off_grid < rates[1] and size > 0.05
    print(
        f"SCPC selected q={q} at decay {off_grid:.6f}: cutoff {cv:.6f}, size {size:.8%}"
    )


def spur_grid_case():
    rng = np.random.default_rng(2971)
    tuning = rng.standard_normal((2, 3000))
    # Recreate the geometry of the residual diagnostic fixture.
    for _ in range(942):
        n = int(rng.integers(6, 15))
        sites = np.cumsum(np.exp(rng.uniform(-7, 2, n)))
        coords = np.column_stack((sites, rng.normal(size=n) * rng.choice([0.2, 2, 20])))
        design = np.column_stack((np.ones(n), rng.normal(size=(n, 2))))
    d = squareform(pdist(coords))
    d /= d.max()
    distances = d[np.tril_indices(n, -1)]
    qx, _ = np.linalg.qr(design)
    m = np.eye(n) - qx @ qx.T
    bm = -0.5 * m @ d @ m
    r = get_r(bm, 2)

    def omega(rho):
        rate = brentq(
            lambda c: np.exp(-c * distances).mean() - rho, 0, 40 / distances.min()
        )
        return r.T @ np.exp(-rate * d) @ r

    null = omega(0.001)
    projected_bm = r.T @ bm @ r
    alternative = (
        null + get_ha_param_i0(null, null, projected_bm, tuning) * projected_bm
    )
    inv_null, inv_alt = np.linalg.inv(null), np.linalg.inv(alternative)
    grid = np.array([omega(rho) for rho in np.linspace(0.0001, 0.03, 30)])
    grid_chol = np.linalg.cholesky(grid)
    limit_chol = np.eye(2)

    def size(chol, cutoff):
        eig = np.linalg.eigvalsh(
            chol.swapaxes(-1, -2) @ (inv_null - cutoff * inv_alt) @ chol
        )
        ratio = np.maximum(eig[..., 1], 0) / np.maximum(
            -eig[..., 0], np.finfo(float).tiny
        )
        return np.where(eig[..., 0] >= 0, 1.0, 2 / np.pi * np.arctan(np.sqrt(ratio)))

    lr_range = np.linalg.eigvals(np.linalg.solve(inv_alt, inv_null)).real
    before = brentq(lambda c: size(grid_chol, c).max() - 0.05, *sorted(lr_range))
    after = brentq(
        lambda c: max(size(grid_chol, c).max(), size(limit_chol, c)) - 0.05,
        *sorted(lr_range),
    )
    assert size(limit_chol, before) > 0.05 and after > before
    print(
        f"SPUR residual I(0) at independence: size {size(limit_chol, before):.8%} → {size(limit_chol, after):.8%}"
    )


if __name__ == "__main__":
    scpc_grid_case()
    spur_grid_case()
