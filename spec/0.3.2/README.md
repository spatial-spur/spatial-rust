# SCPC and SPUR methods specification 0.3.2

The statistical targets are the [conditional OLS interval](../0.1.0/scpc-conditional-ols.md), the [transposed conditional IV weights](../0.2.0/scpc-iv.md) and the [SPUR half-life confidence set](../0.1.0/spur-half-life.md). This record specifies numerical and package conventions for producing those targets. It extends [version 0.3.1](../0.3.1/README.md).

## Selected design and fixed effects

The fitted model supplies the estimation sample, selected coefficient columns, observation-level scores, bread and actual intercept convention. For clustered inference each active observation must have a group ID, and every coordinate must be identical within a group. Sum the scores for each spatial site before applying the bread. Conditional IV uses the structural regressors and full instrument design to define the projector. Its score direction uses the fitted regressors in the estimating equations. Apply the transposed IV residual operator to spatial weights in the order specified in version 0.2.0.

After absorbing fixed effects, divide each nonzero column of the design and instrument matrices by its maximum absolute entry. For a selected OLS design, require as many singular values above $10^{-10}\sigma_{\max}$ as selected columns. For the instrument matrix, keep the left singular directions above the same relative threshold. Redundant instruments are admissible. For IV identification, retain the [projected-design threshold with its absolute floor](../0.3.1/README.md). These decisions apply after the host estimator has selected its model.

For a single absorbed effect, compare the largest absolute within residual in each group and column with half the binary64 spacing at that group's largest absolute input. Zero the group’s residuals only when they fall strictly below this threshold. Use the finite binade spacing at the largest finite double, and the subnormal spacing for subnormal inputs and zero. Compute the comparison in scaled coordinates so half a subnormal spacing is representable during the test. The rule retains a balanced one-spacing contrast. For multiple absorbed effects, use the corresponding column-wide threshold and check the final group moments. Apply absorption to each spatial right-hand side before removing the within-design component. The screening rule is for model columns.

## Conditional weights and rejection probabilities

Multiplying every conditional weight for one coefficient by the same positive constant leaves its Gaussian rejection event invariant. Use relative projection-degeneracy checks and scale the weights together before forming their covariance. Finite changes of units should preserve the fitted subspace, cutoffs and p-values for representable inputs and outputs.

For $G\sim N(0,\Omega)$ the rejection event is $G_0^2 > c^2\sum_{j=1}^q G_j^2$ when $c$ is the cutoff in the unscaled Gaussian ratio. Equality contributes zero, including a degenerate distribution concentrated on the boundary. The covariance must be symmetric and positive semidefinite up to floating-point roundoff. With $q=1$, unit variances and $c=1$, exchangeability gives probability $1/2$ for every $|\rho|<1$. At $\rho=1$ and $\rho=-1$ the probability is zero. A matrix with zero diagonal and unit off-diagonal entries is indefinite. These cases provide provider-independent checks.

## Conditional Nyström calibration

Use one sampled pair chain for the unconditional covariance grid. Starting at the generator state immediately after that chain, draw one second pair chain and use it for **every** conditional coefficient of that fit. Use the same chain for every decay value in its grid. For a fixed model, site ordering and seed, reordering or subsetting the reported coefficients must preserve the inference for each retained coefficient. The first conditional coefficient retains the pair sample that earlier implementations used. This convention affects Monte Carlo calibration in the approximate method.

## SPUR pipeline dispatch

Use one fixed Gaussian draw matrix for all four diagnostics. In automatic mode a tied eigenspace in a native diagnostic sends the four diagnostics and the LBM transformation through the language fallback. Reuse the original draws and report the backend for each stage. The two SCPC fits retain their own dispatch. Preserve the [geometry and simulation conventions](../0.3.0/README.md) and the reported scientific quantities across equivalent geometry storage layouts.
