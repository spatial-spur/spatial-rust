# Spatial SPUR methods version 0.4.0

The shared Rust library owns the numerical methods. Python, R, and Stata own
model fitting, fitted-sample extraction, row identity, and result presentation.
The NumPy/SciPy reference constructs small problems with dense projections and
independent probability integration. The versioned [source register](../0.1.0/sources.md)
identifies the papers and original routines. Historical files record the
conventions of their own versions.

## Support matrix

| Fit or workflow | Python | R | Stata |
| --- | --- | --- | --- |
| Unweighted OLS | Statsmodels OLS and PyFixest feols | lm and fixest feols | regress with the original estimation sample |
| Ordinary unweighted 2SLS | PyFixest Feiv with stored within X and Z | fixest feols with retained structural and first-stage data | ivregress 2sls without absorb |
| Absorbed categorical effects | PyFixest one-way and multiway | fixest one-way and multiway | Unconditional areg slopes with one plain absorb variable |
| Weighted fits and offsets | Unsupported | Unsupported | Weighted fits and unsupported estimation commands are rejected |
| Site clustering | One explicit site label per observation | One explicit site label per observation | cluster() or a single fitted cluster with constant site coordinates |
| SPUR residual tests and pipeline | OLS designs with a constant span | OLS designs with a constant span | Numeric OLS regressors with an intercept |

Conditional SCPC calibrates Gaussian ratio probabilities using the actual fitted
nuisance structure. Unconditional SCPC requires an identified coefficient and a
valid score representation under the relevant weak-dependence assumptions.
Reported intervals are marginal. Conventional 2SLS inference assumes strong
identification, and numerical rank provides only a check on estimability.
Ordinary t intervals and SPUR-based transformation choices carry no general
weak-identification or post-selection guarantee.

## Fitted regression and spatial sites

Keep all estimation observations in fitted order. Validate a mapping from those
rows to supplied coordinate and cluster data before computing weights. Adapters
pass coefficient identities, estimates, actual structural X and instruments Z,
per-observation scores, and bread in coefficient order. Bread is the inverse
total regressor cross product for OLS and the inverse total projected regressor
cross product for 2SLS. It contains no residual variance or degrees-of-freedom
factor. Column aliasing is checked before either conditional or unconditional
calibration. Fixest requires stored within-transformed matrices. Stata validates
the active estimation metadata and materializes the current fitted columns
with the native factor-variable evaluator.

An incidence matrix A maps fitted observations to declared spatial sites. A
site can have repeated fitted observations and one coordinate vector. Without a
cluster label, A is the identity. Different declared sites can have the
same coordinates. They remain separate regression sites and share exact
covariance at the infinite-decay limit.

## SCPC

Let d be the inter-site distance matrix and K(c) have entries exp(-c d).
Calibrate c0 so the mean of its off-diagonal entries equals avc. Let f be
the fraction of co-located pairs. A request needs avc > f. The finite grid has
a 1.2 ratio and includes the endpoint for mean correlation f + (1-f) × 10⁻⁵.
The family also includes K(∞) with entries 1 at equal coordinates and 0
elsewhere. Its projected covariance is the cross product of summed basis rows
at each distinct coordinate.

The spatial working basis consists of the normalized constant W0 and R from
the positive eigenspaces of centered K(c0). Automatic selection evaluates
complete eigenvalue blocks with gaps greater than max(1.6 × 10⁻¹¹,
1.6 × 10⁻¹¹ times the leading magnitude). Positive rank uses a relative
threshold of 10⁻¹². The initial dimension limit is 10, 20, 60, or 120 at avc
cutoffs 0.05, 0.01, and 0.005. Extend the search at a split positive block
or when its largest candidate is selected and more directions are available.
For dimension q, minimize the location-only 95% cutoff times
Γ((q+1)/2)/(√q Γ(q/2)), breaking exact ties toward the smaller q.

For a coefficient with bread-adjusted influence direction a and fitted
residual u, form the site influence h = Aᵀ(a ⊙ u). Its SCPC standard error is
√(n_sites/q) times the norm of Rᵀh. With an absorbed effect, the direction
uses the fitted within-transformed design. For 2SLS the score direction
uses P_Z X. For conditional calibration let s be a normalized a within
each declared site, with the singleton convention s = +1 at zero. The
numerator weight is W0 multiplied by Aᵀ(a ⊙ s). The denominator weights
are Aᵀ diag(s) Mᵀ M_FE diag(a) A R. M_FE applies the fitted fixed-effect
residualizer. OLS uses the orthogonal X residual projection Mᵀ.
2SLS uses Mᵀ = I − (P_Z X)(XᵀP_Z X)⁻¹Xᵀ. For multiway effects,
the host residualizes blocks of spatial directions using the fitted groups.

Obtain size and p-values from Gaussian ratio tails over both the projected
location family and the coefficient-specific family. Unconditional inference
uses the location family. The calculation includes K(∞) for every candidate,
coefficient, critical value, and p-value. The default critical-value tolerance
is 0.001 in t units. Finite search over the interior allows numerical
calibration error and provides no bound over a continuous nuisance parameter.

Exact SCPC uses full geometry. Auto requests Nyström at 4,500 sites. Nyström
samples landmarks in lexicographic coordinate order with a recorded seed.
Euclidean coordinates are centered at the bounding midpoint and divided by
the largest absolute centered coordinate. Geographic points are canonicalized.
Embed complete positive landmark blocks by applying the covariance to
landmark eigenvectors divided by the square roots of their eigenvalues. Use a
common eigenvalue scale within each numerical tie. Use one Frobenius scale per embedded sample to preserve
rotation invariance and take the leading pooled subspace. The complete
landmark sample gives the exact working subspace. The finite covariance grid
is projected by tiles over all sites, including K(∞), so random pair sampling
does not enter calibration. Default landmark settings are 1,000 sites and
20 samples. Review approximate results with exact geometry or different
landmark seeds when scientific conclusions are sensitive to concentrated
leverage. Weight rows are restored to reported site order.

## Distances, spatial diagnostics, and half-life

Euclidean distances use the supplied coordinate units. Geographic distance
is the great-circle arc on a sphere with radius 6,371,008.8 m. Compute
differences in degrees before converting angles. Normalize equivalent
longitudes and use one longitude at an exact pole. SCPC geographic decay
uses fractions of the Earth’s circumference. Geographic ISO radii and
unnormalized half-life results use meters.

SPUR diagnostics use distances divided by the sample diameter. The centered
Lévy–Brownian covariance is the projection of −d/2. The residual version
projects it off the full OLS design, which must span a constant. A fixed q
must fit within positive rank and end at a complete eigenvalue block. Tied
blocks are oriented from their projector with coordinate-ordered anchors.
Use q in the validated residual dimension and interpret results based on the
spatial benchmark, the design and disturbance assumptions, and the number of
spatial contrasts.

I(0) calibrates the angular Gaussian likelihood ratio over a finite grid from
the attainable co-location floor to average correlation 0.03. The alternative
adds a fitted nonnegative Brownian component to a covariance at correlation
0.001 or the co-location limit. I(1) uses the projected Brownian null and
a fitted exponential alternative. Both tune to approximately 50% power on
a separate Gaussian block and report attained power when the target is
unreachable. Use the inference block for p = (1 + number at least as extreme)
 / (1 + draws) and corresponding 1%, 5%, and 10% order-statistic cutoffs.
A size below the minimum attainable p-value has an infinite cutoff.

Half-life uses the angular likelihood
|Σ|⁻¹ᐟ²(zᵀΣ⁻¹z)⁻ᑫᐟ² in the leading Brownian directions. The alternative mixes
exponential covariances with half-life uniform from zero to one sample
diameter. It uses nested panel quadrature and stable log accumulation.
The null includes the exact co-location limit at zero half-life and the
Brownian limit at infinite half-life. Invert finite-simulation tests over a
grid in h/(1+h), with close-pair probes and boundary refinement. Save every
detected confidence component, including gaps and unbounded support. Finite
probes can miss narrow components or gaps. The implementation accepts q up
to 100 for half-life. That limit defines the supported quadrature workload.

## Transformations and execution

LBM-GLS applies the centered Brownian inverse square root to the fitted
variables. The positive whitening spectrum retains eigenvalues above
10⁻¹⁰ at unit diameter. The combined workflow reuses the full decomposition
for its levels diagnostics and whitening. Standalone diagnostics can use
a checked partial eigensolve. Covariance projections are tiled. LBM
eigensolves and exact SCPC still require dense workspaces and may need
quadratic distance work.

Nearest-neighbor differencing averages all other observations at locations
at the smallest positive distance. Numerical ties use 16 binary64 epsilons
times the nearest distance plus a bound for Euclidean input-coordinate
rounding. Exact co-location remains separate. ISO includes the other
observations at distances strictly less than the supplied radius. It returns
zero for an isolated observation. Cluster transformation demeans within
the supplied clusters. NN and ISO use a metric index and centered response
differences. A large radius can still require queries and accumulation over
a dense neighborhood graph.

Execution uses the shared Rust engine with Accelerate on Apple targets and
faer on other targets. Provider and thread choices do not change the method.
Seeded draws are local to a computation. Supplied draws make cross-host
comparisons independent of each language’s RNG. Errors in identification,
calibration, and eigenspace eligibility propagate to the caller. Stata
writes requested outputs only after successful computation.
