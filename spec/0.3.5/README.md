# SCPC and SPUR methods specification 0.3.5

This record extends the implementation conventions in [0.3.4](../0.3.4/README.md). The paper-level targets and author replication profiles are recorded independently in [0.1.0](../0.1.0/conventions.md).

## SCPC

Project the centered coefficient influence values to estimate the variance. Use the fitted coefficient divided by the resulting standard error for the t-statistic against zero. For example, adding a constant to an intercept-only outcome with representable residual variation preserves the standard error.

An exact-zero direction in a singleton cluster uses the unclustered sign convention of $+1$. This preserves conditional weights when each observation is assigned its own cluster. An all-zero direction in a cluster with multiple observations retains the zero vector convention.

Keep the geometric decay candidates and include the finite upper endpoint when the geometric grid ends below it. Both full-matrix and Nyström calculations evaluate a finite nuisance grid. The exponential correlation benchmark and the associated asymptotic assumptions govern coverage. The intervals are marginal for individual coefficients.

The linear conditional calculation requires a linear estimating model. In R that includes `lm`, `fixest::feols` including its linear IV fits, and Gaussian identity-link GLMs supported by the fallback.

Ordinary IV t inference requires strong identification and the relevant score-process assumptions. The conditional 2SLS result takes the larger of the conditional and unconditional cutoffs. Its conditional size argument requires the specified Gaussian benchmark given the full realized design, including the endogenous regressors. With one endogenous regressor $p$ and one excluded instrument $z$, an Anderson–Rubin test of $\beta=b_0$ tests the coefficient of $z$ in the regression of $y-b_0p$ on $z$ and the exogenous controls. Invert those scalar tests for a confidence set. [Müller and Watson (2023), §5.3](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=15)

## SPUR

Residual diagnostics require the nuisance design to span a constant. A full set of group indicators can supply that direction. For $M_X\mathbf{1}=0$, evaluating the projected exponential proxy as $M_X\operatorname{expm1}(-cD)M_X$ removes the canceled all-ones term before multiplication. The decay at average correlation 0.999 sets the basis proxy. When $M_X\mathbf{1}\ne0$, the proxy has an order-one term $(M_X\mathbf{1})(M_X\mathbf{1})^\top$. Extending the residual diagnostics to that case requires defining an anchored model and its calibration. [Müller and Watson (2024), §§2.5 and 4.5](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf)

The combined formula pipeline requires an explicit intercept. LBM GLS centers the transformed columns. The intercept in a transformed fit reflects that normalization. In manual transformed IV, the analyst needs exclusion and relevance for the transformed moments, since the transformation can mix instruments and errors across sites.

Orient each simulated basis vector by the first entry whose magnitude is within a relative $10^{-12}$ of its largest magnitude. This stabilizes the sign of nearly equal loadings across hosts.

Decay calibration requires a target strictly above the co-location correlation floor and below one. Verify the bracket before bisection. The ordinary, bracketed calculation retains the legacy width tolerance. Continue the search and refine to relative accuracy when the root lies beyond the legacy upper trial. The nearest-neighbor fallback compares raw positive distances before averaging any ties.

Half-life reports the hull of accepted points on the [Table 1 replication grid](../0.1.0/conventions.md#author-table-1-replication). The last tested point is finite at 100 sample diameters. An accepted last point is displayed as an open upper endpoint under this reporting convention. The published continuous inversion has its own [independent reference calculation](../0.1.0/reference/README.md#spur-half-life-continuous-inversion).

The native bindings use SCPC model API 7 and SPUR API 4.
