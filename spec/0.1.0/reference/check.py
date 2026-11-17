"""Check the analytic fixtures with dense linear algebra and independent integration."""

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import t as student_t


def centered_basis(covariance, q):
    n = covariance.shape[0]
    center = np.eye(n) - np.ones((n, n)) / n
    eigenvalues, eigenvectors = np.linalg.eigh(center @ covariance @ center)
    return eigenvalues[-q:][::-1], eigenvectors[:, -q:][:, ::-1]


def conditional_scpc(case):
    x = np.asarray(case["regressor"], dtype=float)
    y = np.asarray(case["outcome"], dtype=float)
    n = len(y)
    q = case["q"]
    design = np.column_stack((np.ones(n), x))
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    residuals = y - design @ coefficients
    residual_projection = np.eye(n) - design @ np.linalg.pinv(design)
    sites = np.asarray(case["sites"], dtype=float)
    distances = np.abs(sites[:, None] - sites)
    off_diagonal = distances[~np.eye(n, dtype=bool)]
    decay = brentq(
        lambda c: np.exp(-c * off_diagonal).mean() - case["max_average_correlation"],
        0,
        2,
        xtol=1e-14,
    )
    kernel = np.exp(-decay * distances)
    _, basis = centered_basis(kernel, q)
    constant = np.ones(n) / math.sqrt(n)
    signs = np.sign(x)[:, None]
    x_column = x[:, None]

    # The proof in README.md locates both continuous suprema at this decay.
    loc = np.column_stack((constant, basis))
    cond = np.column_stack(
        (
            abs(x) * constant,
            signs * (residual_projection @ (x_column * basis)),
        )
    )
    historical = np.column_stack(
        (
            cond[:, 0],
            x_column * (residual_projection @ (signs * basis)),
        )
    )
    loc_cov, cond_cov, historic_cov = (
        weights.T @ kernel @ weights for weights in (loc, cond, historical)
    )
    for covariance in (loc_cov, cond_cov, historic_cov):
        np.testing.assert_allclose(covariance[0, 1], 0, atol=1e-12, rtol=0)

    se = math.sqrt(n) * np.linalg.norm(basis.T @ (x * residuals))
    se /= float(x @ x) * math.sqrt(q)
    beta = float(coefficients[1])
    tau = (beta - case["null_coefficient"]) / se
    loc_scale, cond_scale, historic_scale = (
        math.sqrt(cov[0, 0] / cov[1, 1]) for cov in (loc_cov, cond_cov, historic_cov)
    )
    factor = 1 / math.tan(math.pi * case["alpha"] / 2)
    critical = max(loc_scale, cond_scale) * factor
    return {
        "decay": decay,
        "coefficient": beta,
        "se": se,
        "t": tau,
        "location_covariance": loc_cov,
        "conditional_covariance": cond_cov,
        "historical_covariance_at_c_min": historic_cov,
        "location_p": 2 * math.atan2(loc_scale, abs(tau)) / math.pi,
        "conditional_p": 2 * math.atan2(cond_scale, abs(tau)) / math.pi,
        "historical_conditional_p_at_c_min": 2
        * math.atan2(historic_scale, abs(tau))
        / math.pi,
        "location_critical_value": loc_scale * factor,
        "conditional_critical_value": cond_scale * factor,
        "interval": [beta - critical * se, beta + critical * se],
    }


def scpc_component_selection(case):
    y = np.asarray(case["outcome"], dtype=float)
    distances = np.asarray(case["distance_matrix"], dtype=float)
    n = len(y)
    rho = case["max_average_correlation"]
    kernel = np.exp(math.log(rho) * distances)
    _, basis = centered_basis(kernel, n - 1)
    # Increasing variance ratio with rho gives the worst case at this endpoint.
    scale = math.sqrt((1 + (n - 1) * rho) / (1 - rho))
    counts = range(1, n)
    critical = [scale * student_t.isf(case["alpha"] / 2, q) for q in counts]
    lengths = [
        math.sqrt(8 / n)
        * cutoff
        * math.gamma((q + 1) / 2)
        / (math.sqrt(q) * math.gamma(q / 2))
        for q, cutoff in zip(counts, critical, strict=True)
    ]
    q = int(np.argmin(lengths)) + 1
    beta = float(y.mean())
    se = np.linalg.norm(basis[:, :q].T @ (y - beta)) / math.sqrt(n * q)
    tau = (beta - case["null_coefficient"]) / se
    width = critical[q - 1] * se
    return {
        "candidate_critical_values": critical,
        "expected_iid_interval_lengths": lengths,
        "selected_q": q,
        "coefficient": beta,
        "se": se,
        "t": tau,
        "p": 2 * student_t.sf(abs(tau) / scale, q),
        "interval": [beta - width, beta + width],
    }


def spur_half_life(case):
    sites = np.asarray(case["sites"], dtype=float)
    distances = np.abs(sites[:, None] - sites)
    diameter = float(distances.max())
    eigenvalues, basis = centered_basis(-distances / 2, case["q"])
    z = basis.T @ case["outcome"]
    first_direction_share = float(z[0] ** 2 / (z @ z))
    c = math.log(2) / diameter
    null_cov = basis.T @ (np.exp(-c * distances) / (2 * c)) @ basis

    def covariance_ratio(half_life):
        rho = 0 if half_life == 0 else math.exp(-math.log(2) / half_life)
        return 3 * (1 + rho) / (3 - rho)

    def density(half_life):
        ratio = covariance_ratio(half_life)
        return math.sqrt(ratio) / (
            first_direction_share + ratio * (1 - first_direction_share)
        )

    # The proof in README.md gives one-sided angular tails for h >= diameter.
    def null_p(half_life, observed):
        ratio = covariance_ratio(half_life)
        return (
            2 * math.atan(abs(observed[0] / observed[1]) / math.sqrt(ratio)) / math.pi
        )

    mixture = quad(density, 0, diameter, epsabs=1e-13, epsrel=1e-13)[0] / diameter
    upper = brentq(lambda h: null_p(h, z) - case["alpha"], diameter, 128, xtol=1e-13)
    unbounded_z = basis.T @ case["unbounded_case_outcome"]
    return {
        "lbm_eigenvalues": eigenvalues,
        "null_covariance_at_diameter": null_cov,
        "mixture_density_without_common_constant": mixture,
        "likelihood_ratio_at_diameter": mixture / density(diameter),
        "null_p_values": [
            null_p(math.inf if h == "infinity" else h, z) for h in case["tested_nulls"]
        ],
        "upper_endpoint": upper,
        "normalized_upper_endpoint": upper / diameter,
        "unbounded_case_p_at_infinity": null_p(math.inf, unbounded_z),
    }


def main():
    path = Path(__file__).with_name("fixtures.json")
    fixtures = json.loads(path.read_text(encoding="utf-8"))
    if fixtures["spec_version"] != path.parents[1].name:
        raise ValueError("fixture version and spec directory differ")
    tolerance = fixtures["absolute_tolerance"]
    for calculate in (conditional_scpc, scpc_component_selection, spur_half_life):
        name = calculate.__name__
        case = fixtures[name]
        actual = calculate(case)
        expected = case["expected"]
        if actual.keys() != expected.keys():
            raise ValueError(f"{name}: reference and fixture fields differ")
        for key, value in actual.items():
            np.testing.assert_allclose(
                value,
                expected[key],
                atol=tolerance,
                rtol=0,
                err_msg=f"{name}.{key}",
            )
        print(f"{name}: {len(expected)} checks passed")


if __name__ == "__main__":
    main()
