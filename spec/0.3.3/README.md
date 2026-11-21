# SCPC and SPUR methods specification 0.3.3

This record extends the numerical and host conventions in [0.3.2](../0.3.2/README.md). The [conditional OLS interval](../0.1.0/scpc-conditional-ols.md), [IV projection](../0.2.0/scpc-iv.md) and [SPUR half-life confidence set](../0.1.0/spur-half-life.md) continue to define the statistical targets.

## Numerical limits

At a zero cutoff the Gaussian rejection probability is one when the numerator variance is positive and zero when the numerator is identically zero. For positive cutoffs use the strict quadratic inequality in 0.3.2. Compare covariance scales before squaring the cutoff so a finite variance ratio stays representable. The diagonal case with numerator variance $10^{200}$, two denominator variances $10^{-140}$ and cutoff $10^{170}$ gives $1-1/\sqrt{2}$.

Scale Euclidean coordinates before constructing the cross-product used by Nyström's canonical rotation. For a fixed sample with isolated principal directions, changing units of all coordinates preserves the canonical order. If the component search reaches its current limit, increase the limit to $\lceil 1.5 q_{\max}\rceil$ before applying the method's upper bound.

## Model and formula handling

Conditional SCPC absorbs the entire fitted nuisance space. For a fitted varying slope this includes the corresponding group-by-variable columns and any fitted group intercepts. The R fallback uses fitted group and slope information. A varying-slope model with unavailable fitted data yields an explicit error.

Evaluate single-series SPUR formulas on the estimation sample. Preserve expression semantics and the fitted design's intercept convention. Half-life likelihood comparisons preserve the strict ordering of the likelihood ratios when terms approach underflow. The host fallbacks use log-space comparisons with the same Gaussian draws.
