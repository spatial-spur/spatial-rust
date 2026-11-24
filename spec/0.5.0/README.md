# Spatial SPUR methods version 0.5.0

The shared Rust library owns the numerical methods. Python, R, and Stata own
model fitting, fitted-sample extraction, row identity, and result presentation.
The NumPy/SciPy reference constructs small problems with dense projections and
independent probability integration. The versioned [source register](../0.1.0/sources.md)
identifies the papers and original routines. Historical files record the
conventions of their own versions.

## Support matrix

| Fit or workflow | Python | R | Stata |
| --- | --- | --- | --- |
| Unweighted OLS | Statsmodels OLS and PyFixest feols | lm and fixest feols | regress with verified estimation sample & normal equations |
| Ordinary unweighted 2SLS | PyFixest Feiv with stored within X and Z | fixest feols with retained structural and first-stage data | ivregress 2sls without absorb |
| Absorbed categorical effects | PyFixest one-way and multiway | fixest one-way and multiway | Unconditional areg slopes with one plain absorb variable |
| Weighted fits and offsets | Unsupported | Unsupported | Weighted fits and unsupported estimation commands are rejected |
| Site clustering | One explicit site label per observation | One explicit site label per observation | cluster() or a single fitted cluster with constant site coordinates |
| SPUR residual tests and pipeline | OLS designs with a constant span | OLS designs with a constant span | Numeric OLS regressors with an intercept |

Conditional SCPC calibrates Gaussian ratio probabilities using the actual fitted
nuisance structure. Unconditional SCPC requires an identified coefficient and a
valid score representation under the relevant weak-dependence assumptions.
Reported intervals are marginal.

**Conditional 2SLS scope:** Valid 2SLS t-statistics and confidence intervals rely
on strong first-stage identification; numerical rank only checks estimability,
not instrument strength. Under endogeneity, the finite-sample conditional
Gaussian argument is heuristic (SCPC protects asymptotic validity by taking the
maximum over unconditional and conditional critical values). For weak
instruments, scalar Anderson–Rubin testing on structural null hypotheses is the
appropriate alternative; multi-parameter joint or weak-IV inference is outside
the default scalar SCPC scope. Similarly, SPUR-selected transformations carry no
post-selection uniformity guarantee.

## Fitted regression and spatial sites

Keep all estimation observations in fitted order. Validate a mapping from those
rows to supplied coordinate and cluster data before computing weights. Do not
accept positional-only external frame matches from non-identifying shared
columns. For R `fixest`, apply all sample-selection indices (`obsRemoved`, split
selections) after `subset`. Stata verifies both residual sum of squares and
regressor/instrument normal equations before accepting an active fit.

An incidence matrix A maps fitted observations to declared spatial sites. A
site can have repeated fitted observations and one coordinate vector. Without a
cluster label, A is the identity. Different declared sites can have the same
coordinates. They remain separate regression sites and share exact covariance at
the infinite-decay limit.

## SCPC

Let d be the inter-site distance matrix and K(c) have entries exp(-c d).
Calibrate c0 so the mean of its off-diagonal entries equals avc. Let f be the
fraction of co-located pairs. A request needs avc > f. The finite upper
correlation target is `f + min((1 - f) * 10⁻⁵, 0.1 * (avc - f))`, guaranteeing
strictly positive grid span (`cmax > c0`). The finite grid uses a **1.05**
geometric ratio (`CGRID_FAC = 1.05`) from `c0` through `cmax` plus `K(∞)`.

The spatial working basis consists of the normalized constant W0 and R from
the positive eigenspaces of centered K(c0). Automatic selection evaluates
complete eigenvalue blocks with gaps greater than max(1.6 × 10⁻¹¹,
1.6 × 10⁻¹¹ times the leading magnitude). Positive rank uses a relative
threshold of 10⁻¹². The initial dimension limit is 10, 20, 60, or 120 at avc
cutoffs 0.05, 0.01, and 0.005. Extend the search at a split positive block or
when its largest candidate is selected and more directions are available (even in
single-sample Nyström). For dimension q, minimize the location-only 95% cutoff
times Γ((q+1)/2)/(√q Γ(q/2)), breaking exact ties toward the smaller q.

Exact SCPC uses full geometry. Auto requests Nyström at 4,500 sites. Nyström
calibrates `[c0, cmax]` over all full-site pairs (streamed in `O(n)` memory),
samples landmarks in canonical lexicographic `(lat, lon)` / Euclidean order
with a recorded seed, skips degenerate all-co-located subsets whenever globals
allow `avc > f`, and re-centers/orthonormalizes pooled embedded directions.

## Distances, spatial diagnostics, and half-life

Euclidean distances use Euclidean coordinates; geographic distance uses great-circle
chordal / arc metrics normalized as fraction of Earth circumference (SCPC) or
meters (unnormalized half-life / ISO). SPUR tests divide distances by the sample
diameter. Both simple and tied eigenvectors of the Lévy–Brownian basis are
oriented using canonical coordinate-ordered anchors so fixed-draw results are
invariant to observation row permutation.

I(0) calibrates the angular Gaussian likelihood ratio over a hybrid log+linear
80-point correlation grid over `(floor, 0.03]` (48 log points near `floor` plus
32 linear points through `0.03`). Both I(0) and I(1) always report the achieved
best tuning power (`best_tuning_power`), and all simulation draw columns must
have finite strictly positive L2 norms (all-zero draws are rejected).

Half-life uses 10-point Gauss–Legendre composite dyadic panel quadrature on
`[0, 1]` and inverts tests over a 256-step grid in `h/(1+h)` augmented with
dense close-pair (`h -> 0`) and large-scale (`h -> ∞`) logarithmic probes, plus
adaptive multi-round midpoint insertion on borderline accepted/rejected intervals
to capture narrow gaps or islands prior to root bisection.

## Transformations and execution

In Python (`patsy`) and R (`model.frame`), `spurtransform()` evaluates formula
terms/expressions (rejecting `offset(...)`) on complete observations and pads
incomplete rows with `NaN` / `NA`. The combined `spur()` pipeline evaluates the
LBM-GLS transformed regression **without an intercept** (`nocons`), reporting
only the spatially whitened slope regressors across Python, R, and Stata.

Native R bindings use only the official public R C-API without `extendr-api`.
