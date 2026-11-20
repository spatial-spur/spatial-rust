# SCPC and SPUR methods specification 0.3.1

This version extends the numerical safeguards in the [0.3.0 implementation conventions](../0.3.0/README.md). The [transposed conditional IV weights](../0.2.0/scpc-iv.md) still apply.

## Numerical identification for conditional 2SLS

After removing any absorbed effects, equilibrate each nonzero regressor column to unit maximum absolute entry. Let $X_s$ be that design and $Q_Z$ an orthonormal basis for the instrument span. Compute the singular values of $Q_Z' X_s$. Require a full column rank at the threshold

$$
\sigma_{\min}(Q_Z' X_s) > 10^{-10} \max\{1, \sigma_{\max}(Q_Z' X_s)\}.
$$

The positive floor detects a projected design made entirely of roundoff, including when there is only one regressor. The numerical check governs construction of the projector. Empirical instrument strength requires its own econometric assessment.

For absorbed effects, preserve representable within-group variation and prevent residual rounding error from turning an absorbed column into a usable direction during equilibration. The R and Python implementations compare the largest within-column residual to the floating-point spacing at the largest absolute input entry. Keep a residual at or above that spacing. Apply the fixed-effect projection to each spatial right-hand side before removing the within-design component, as required by the projection order. An absorbed instrument can be discarded if the remaining instrument span identifies the design.

Regression checks cover a single nearly orthogonal instrument, a design with all projected directions near zero, weakly connected fixed effects and translated regressors.
