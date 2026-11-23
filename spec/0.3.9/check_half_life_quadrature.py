"""Check the half-life mixture quadrature on a regular line, independently of SPUR."""

import json

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp, roots_legendre


def log_mixtures(order):
    n = 128
    sites = np.linspace(0, 1, n)
    distances = np.abs(sites[:, None] - sites)
    dimensions = (4, 15, 40, 100)
    modes = np.arange(1, max(dimensions) + 1)
    # The nonconstant eigenvectors of the centered Brownian covariance on a line.
    basis = np.sqrt(2 / n) * np.cos(np.pi * (np.arange(n)[:, None] + 0.5) * modes / n)
    directions = {}
    for q in dimensions:
        vectors = np.column_stack(
            (np.eye(q)[:, 0], np.eye(q)[:, -1], (-1.0) ** np.arange(q))
        )
        directions[q] = vectors / np.linalg.norm(vectors, axis=0)

    first_exponent = int(np.floor(np.log2(min(2.0**-9, 1 / (32 * (n - 1))))))
    edges = np.r_[0, 2.0 ** np.arange(first_exponent, 1)]
    nodes, weights = roots_legendre(order)
    log_terms = {q: [] for q in dimensions}
    for left, right in zip(edges[:-1], edges[1:]):
        half_width = (right - left) / 2
        for node, weight in zip(nodes, weights):
            h = left + half_width * (node + 1)
            covariance = basis.T @ (np.exp(-np.log(2) * distances / h) @ basis)
            lower = np.linalg.cholesky((covariance + covariance.T) / 2)
            for q in dimensions:
                factor = lower[:q, :q]
                whitened = solve_triangular(factor, directions[q], lower=True)
                log_density = -np.log(np.diag(factor)).sum() - (q / 2) * np.log(
                    (whitened**2).sum(axis=0)
                )
                log_terms[q].append(log_density + np.log(weight * half_width))
    return {q: logsumexp(terms, axis=0) for q, terms in log_terms.items()}


if __name__ == "__main__":
    production, check, reference = (log_mixtures(order) for order in (5, 20, 32))
    results = {}
    for q in production:
        errors = np.abs(np.expm1(production[q] - reference[q]))
        convergence = float(np.max(np.abs(np.expm1(check[q] - reference[q]))))
        assert max(errors) < 1e-5 and convergence < 1e-10, (q, errors, convergence)
        results[q] = {
            "relative_error_5_to_32": errors.tolist(),
            "max_relative_difference_20_to_32": convergence,
        }
    print(json.dumps(results, indent=2))
