# SCPC and SPUR methods 0.3.9

This revision inherits the [0.3.8 record](../0.3.8/README.md).

## Statistical definitions

For ISO differencing with radius \(b>0\), let \(N_i=\{j:j\ne i,\ d_{ij}<b\}\). The transformed value is \(y_i-|N_i|^{-1}\sum_{j\in N_i}y_j\) when \(N_i\) is nonempty. Distinct observations at the same site belong to each other's neighborhoods. A pair exactly \(b\) apart is excluded. The package assigns zero to an isolated row. Nearest-neighbor differencing continues to average ties at the smallest positive distance. [Müller and Watson, “Spatial Unit Roots”, spatial transformations](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf#page=21)

Write \(M=I-\mathbf1\mathbf1'/n\) and \(K=M(-D/2)M\) using distances normalized by the sample diameter. The LBM-GLS transformation on a column \(y\) is \(M K_{\mathrm{ret}}^{-1/2} M y\), where the inverse square root acts on the retained positive eigenspace. The numerical implementation centers each input before whitening and recenters each output. Adding a representable constant to a column preserves the contrasts. The previous eigenvalue retention threshold still applies.

For conditional 2SLS, let \(\widetilde X,\widetilde Z\) denote the structural regressors and instruments after absorbing the fitted fixed effects. The fitted first-stage directions are \(P_{\widetilde Z}\widetilde X\). The directions use the retained regression sample and fitted transformation state, including spline knots. An endogenous regressor can instrument itself. Fixed-effect absorption and the conditional SCPC projection retain their existing order.

Residual SPUR requires the design to span the constant. A nearly constant regressor on its own can leave a genuine constant component outside the fitted space.

## Numerical conventions and user diagnostics

- Test the span of the unit constant after scaling each design column by its largest absolute entry. A projected residual norm up to \(8\epsilon\) is roundoff. Otherwise use the backward-error allowance \(64\epsilon(1+\|A\|_2\|A^+(\mathbf1/\sqrt n)\|_2)\) for the scaled design \(A\). This accounts for cancelling coefficients. Retain the existing response projection correction on accepted designs above the roundoff threshold.
- SPUR requires integer \(2\le q<n\), \(q\le n-\operatorname{rank}(X)\) for residual diagnostics, sufficient positive covariance rank, and a gap at the cutoff. A complete retained tied block is admissible and reported. A split block is ambiguous. Rank, gap and tie tolerances are inherited from 0.3.8.
- Give a single advisory per call when \(q\ge40\) or all numerically positive directions are retained. Higher spatial frequencies can increase sensitivity to local noise. For half-life inference, \(2\le q\le5\) warrants a low-information note. The fixed five-node-per-panel mixture quadrature is currently supported through \(q=100\). The ceiling describes the numerical implementation. It does not bound the theoretical statistic or the dimension of I(0), I(1), or SCPC.
- Tune the alternative on the independent tuning draws. An exact 50% estimate ends the search, and a bracket found while decreasing the parameter can be reused. If the high-side search fails to reach the 50% target, retain the tested alternative with highest estimated power and report that estimate. For I(0), stop when the trace-normalized mixture is within \(8\epsilon\) entrywise of the limiting Brownian shape. I(1) already stops when every positive-distance correlation rounds to zero. A low-\(q\) diagnostic that attains the target needs no tuning warning.
- With \(B\) inference draws, the smallest rank p-value is \(1/(B+1)\). List the unattainable nominal rejection levels among 1%, 5% and 10%. When \(B/(B+1)\) is below the requested confidence level, every half-life must be accepted. Equality still permits rejection.
- Co-location is evaluated on the effective analysis sites. I(0) requires the co-located-pair fraction \(f\le0.03\). I(1) and half-life inference can use geometries above that fraction if the projected experiment is informative. SCPC requires \(\mathrm{avc}>f\). Approximate SCPC uses the pooled landmark-pair fraction for decay calibration. The packages report the basis for the fraction in an eligibility error.

On 128 equally spaced sites, comparing five and 32 nodes per panel for three response directions gives maximum relative mixture-density differences of \(7.4\times10^{-9}\) at \(q=4\), \(2.0\times10^{-8}\) at \(q=15\), \(9.5\times10^{-8}\) at \(q=40\), and \(3.1\times10^{-6}\) at \(q=100\), rounded up. Twenty and 32 nodes agree to within \(1.2\times10^{-13}\) on these cases. Accuracy can vary with the geometry and response direction. [Reproduce the comparison](check_half_life_quadrature.py)

Native adapters use SCPC API 7 and SPUR API 6. SPUR returns the retained-tie flag, whether all positive directions were used, and the best independent tuning power when the high-side target was unattainable.
