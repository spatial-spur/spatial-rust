# SCPC and SPUR methods specification 0.3.0

This version amends the [implementation conventions in version 0.2.0](../0.2.0/implementation-conventions.md).
The [conditional OLS projection order](../0.1.0/scpc-conditional-ols.md#required-projection-order)
and the [transposed IV weights](../0.2.0/scpc-iv.md) continue to govern every backend.

## SCPC

- For an unclustered coefficient, write its signed score direction as
  \(S_x=\operatorname{diag}(s_i)\), where \(s_i=-1\) for \(x_i<0\) and
  \(s_i=+1\) for \(x_i\geq0\). Apply this convention to both IEEE signs of exact zero.
  This reproduces the author's unclustered MATLAB convention for the diagonal
  sign matrix. A tiny nonzero score keeps its own sign. Clustered SCPC
  continues to use the within-cluster direction normalization.
- Form IV projection factors from a rank-revealing decomposition of the
  projected design after scaling nonzero columns. The identifying condition
  concerns the rank of the projected design. The decomposition should avoid
  normal equations, which square its condition number.
- Reject fits with regression weights until a weighted SCPC profile is defined.
  The rejection applies across backends and includes explicit constant weights.
- The sampled Nyström covariance estimator can yield an indefinite matrix.
  Validate each covariance on the block used for inference, including
  conditional blocks. Allow negative roundoff at dimension-scaled machine
  precision relative to the largest absolute entry of the block. Material
  indefiniteness terminates the operation with a dedicated covariance error.
  Selecting an exact calculation explicitly, when feasible, uses the exact
  covariance integral. Tied Nyström eigenspaces remain a distinct outcome.

## SPUR

- SPUR pipeline diagnostics and levels and LBM-GLS inference use the levels
  regression's estimation sample. The transformation acts on its evaluated
  response and numeric design columns. The transformed regression retains
  the levels fit's intercept convention.
- Generate one common Gaussian matrix for the diagnostics. Reuse it if automatic
  dispatch retries any diagnosis in the host language. Record the backend used
  for the diagnostics, transformation, and SCPC fits when reporting the full
  pipeline.
- When orienting isolated eigenvectors, consider entries within
  \(64\epsilon_{\mathrm{mach}}\) times the vector's largest absolute entry
  tied in magnitude. Select the earliest such entry and make it positive.
  This stabilizes the sign convention for symmetric geometries. The rule
  concerns orientation of individual eigenvectors. A repeated eigenvalue
  still leaves rotations within its eigenspace undetermined.
- Distinguish an empty half-life confidence set from a set that accepts the
  open-upper grid endpoint. A finite set reports both grid endpoints in the
  requested units.
- Let \(m\) be the number of simulated statistics strictly greater than the
  observed statistic out of \(B\) draws. Report \(p=m/B\) and retain a grid
  point exactly when \((B-m)/B < \ell\), where \(\ell\) is the requested
  confidence level. Compute this complement from the counts so rounding
  \(1-\ell\) cannot include a grid point at the nominal equality boundary.

The Python adapter requires SCPC API 6 and SPUR API 3. The corresponding
Stata plugin probes are `ping_v6` and `ping_spur_v3`.
