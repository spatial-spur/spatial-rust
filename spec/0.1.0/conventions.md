# Conventions and decisions

Methods specification 0.1.0. The rows below identify named sources of finite behavior. A recorded author or release profile supports reproduction of that source. Product choices listed as open require separate approval before changing defaults. The mathematical targets are in [Conditional OLS SCPC](scpc-conditional-ols.md) and [SPUR half-life](spur-half-life.md).

## SCPC

| Choice | Author MATLAB in the SPUR archive | Author Stata and pinned package behavior |
| --- | --- | --- |
| Conditional denominator columns | $S_xM_VD_xB_q$, the [required order](scpc-conditional-ols.md#required-projection-order). `get_cv_eta.m` multiplies by `x`, residualizes on the full design, and then multiplies by sign. | The pinned Stata, Python and R paths and Milestone 1 core use $D_xM_VS_xB_q$. This historical defect can change both the critical value and null rejection probability. |
| Sign at zero | `x >= 0` versus `x < 0` in MATLAB, so the sign multiplier is $+1$ at zero. | Ordinary `sign` gives zero. With $x_i=0$ the numerator contribution from that row is zero, but the error can still enter the denominator through $M_V$ under a model allowing it. |
| Correlation bound and finite decay grid | In exact mode, solve for $c_{\min}$ and for the weak-correlation endpoint with target `1e-5`. With multiplier `1.2`, use `m=max(2,ceil(log(cmax/cmin)/log(1.2)))` scales $c_{\min}1.2^j$, $j=0,\ldots,m-1$. | The pinned Stata and Python/R exact paths use the same geometric grid convention. With repeated locations a fixed target of `1e-5` can lie below the attainable floor. The release solver can stop at its scale cap with a missed target. |
| Number of components | Starts the candidate search at 10, 20, 60 or 120 according to the correlation bound. Capped at $n-1$ and expanded by a rounded factor of 1.5 when the winner reaches the current cap. Minimizes the location expected-length criterion at the requested main confidence level, holding the selected $q$ for additional levels. | Stata, Python and R select $q$ at 95% for ordinary inference. The finite cap expansion and rounding belong to the profile. |
| Finite critical values | The archived `findcv.m` starts at a Student-$t_q$ cutoff, searches in units of $k/\sqrt q$, bisects at width `0.001/sqrt(q)`, and reports the last midpoint in ordinary t units. | Author Stata returns the upper bracket endpoint in t units, as do the pinned Python and R implementations in regular cases. The common final conditional cutoff is the maximum of the location and conditional cutoffs. |
| Public p-value | The archived main MATLAB entry point stores `pvalue_scpc`, based on the location family, alongside the conditional interval. | In conditional mode the pinned Stata, Python, and R paths report the maximum of the two family p-values. |

The source locations for these rows are `cscpc/cscpc.m`, `cscpc/get_cv_eta.m`, `cscpc/form_omega_mat.m`, `cscpc/findqcv.m`, and `cscpc/findcv.m` in the [author archive](sources.md#author-software-and-replication), plus author Stata `src/scpc.ado` and the [package revisions](sources.md#package-comparison-revisions). The [four-site reference](reference/README.md#conditional-scpc-with-a-fixed-component-count) independently distinguishes the two conditional weight orders.

The standalone JBES replication materials independently confirm $S_xM_VD_xB_q$ in the MATLAB simulation helper and in the Fortran calculation for nonzero regressors. Their simulation drivers set their own calibration parameters. For example, `Matlab/mc_unitline_model_3_dm.m` uses a weak-correlation endpoint target of `1e-4`. [JBES replication sources](sources.md#author-software-and-replication)

### Open SCPC product choices

- Choose the zero-sign, tied-prefix, nonidentifiable design, and degenerate Gaussian policies. Where $q$ splits a tied eigenspace, record the accepted basis rule or an explicit ambiguity outcome.
- Specify how the software approximates the continuum of decays and verifies the attained rejection bound. Decide which historical grid to provide. An attainable floor-adjusted weak-correlation endpoint is a candidate for a finite grid, while the theoretical supremum is over all $c\ge c_{\min}$.
- Keep numerical tolerances for size, calibration, and probability accuracy in the verification contract. A solver's stopping rule is separate from the statistical target. In particular, a cap in coordinate-dependent units needs its own failure check.

## SPUR half-life

### Author Table 1 replication

- Divide distances by their sample diameter. Express the null and alternative grids as fractions of that diameter. Set $q=15$, the draw count to 100,000, and initialize the Matlab RNG with `rng(197347)` in the Table 1 driver. The seed is for that replication.
- For basis construction, solve for the decay with average pairwise correlation 0.999 and take the leading eigenvectors of its centered exponential covariance. This approximates the paper's centered Lévy–Brownian eigenspace. The remaining covariance calculations use the chosen basis.
- Use 50 equally weighted alternatives `linspace(0.001, 1, 50)`. The null candidates are `linspace(0.001, 1, 100)`, `linspace(1.01, 3, 30)`, and `100`. Each `linspace` includes both endpoints. The candidate `100` is evaluated at a finite normalized half-life.
- Reuse one $q\times B$ matrix of standard normal draws `emat` across nulls and recolor it by each null covariance factor before evaluating the score. Average alternative angular **densities**. Estimate each p-value by the strict exceedance fraction $\#\{L_{\mathrm{sim}}>L_{\mathrm{obs}}\}/B$. The driver keeps candidates with $\widehat p>0.05$, reports the hull, and formats an accepted upper candidate of 100 as infinity. Null values in the gap from 3 to 100 are absent from this reporting profile.

The source locations are the `Matlab/Chetty_Regressions_Table_1_ComputeResults.m` driver, `Matlab/functions/c_ci.m` and the Table 1 formatter. The first three rows of the archived saved normalized intervals are approximately $[0.0918181818,100]$ for absolute mobility, $[0.0413636364,100]$ for fraction Black, and $[0.001,0.2835454545]$ for racial segregation. The formatted table displays infinity in the first two rows. These values document author output from the archive. The archive programs were inspected, and this spec's mathematical test cases were calculated independently. [Archive source paths and digest](sources.md#author-software-and-replication)

### Pinned package behavior

The pinned Python, R and Mata half-life calculations use the same 131/50 grid structure, exponential basis proxy, strict exceedance count, and strict grid acceptance. They generate draws in their respective hosts. For cross-language finite-profile checks use a supplied draw matrix and record basis orientation and layout, since integer seeds select host-specific draws. Python represents the open upper endpoint as infinity, R as `NA`, and Stata as numeric missing. The pinned R and Stata public wrappers apply an Earth-circumference conversion even to Euclidean diameters. Pinned Python preserves Euclidean distance units. This inconsistency belongs to release compatibility. [Package revisions and code locations](sources.md#package-comparison-revisions)

### Open SPUR product choices

- Decide which profile is shipped for half-life inference. The basis proxy and finite grid reproduce the author calculation. The exact LBM basis and continuous inversion implement the published population target.
- Specify whether the result carries the accepted set or only its hull, and how to distinguish empty, bounded, and open-upper sets from an undefined or uninformative experiment. A finite empirical strict p-value rule may reject at equality or with all ratios tied while the population quantile rule retains the null.
- Fix public distance semantics. Dimensionless half-life is a fraction of diameter. Finite absolute half-life uses the units of $D$, including metres for a geographic profile that explicitly defines distances in metres.
- Treat the draw count, common draw policy, RNG and any finite-Monte-Carlo correction as reproducibility choices. Record any changes to the 50-point integration or endpoint sentinel with expected effects and reference tests.
