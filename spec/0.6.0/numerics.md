# Numerical profile 0.6.0

This chapter is normative under the [document map](README.md#1-authority-and-document-map).
Its profile ID is **spatial-spur-finite-mesh-0.6.0**. **[D]**
It defines a practical finite numerical calculation and its diagnostics.
It does not certify a continuous nuisance supremum or complete half-life
inversion. An alternative numerical profile needs its own identifier,
validation, and explicit result metadata.

## 1. Numerical targets and accuracy layers

Keep the following layers distinct:

| Layer | Required target | What it does not establish |
| --- | --- | --- |
| Geometry and linear algebra | Declared distances, projections and subspaces within verified floating-point accuracy | Scientific model assumptions |
| Fixed-covariance probability | Gaussian tail at the given projected covariance, or the declared finite-simulation rank | Nuisance supremum between covariance points |
| SCPC critical-value solve | Upper endpoint of a valid bracket over the same evaluated family as the p-value | A bound on omitted nuisance values |
| I(0) nuisance calibration | Maximum rank p-value/cutoff over the declared finite mesh | Continuum coverage from grid spacing |
| Half-life point test | Correct fixed-null simulation with the declared mixture statistic | Completeness of confidence-set search |
| Half-life search | The bounded recipe and all evaluated points/brackets | Absence of unseen gaps or islands |

The word exact in an SCPC basis option means full geometry/eigenspace
construction. It MUST NOT be presented as exact continuous calibration.
Uniform and adaptive meshes require the same qualification unless a genuine
continuous error bound is separately established.

## 2. Arithmetic, projection, and factorization

Use IEEE binary64 or a demonstrably equivalent/more accurate implementation.
Scale norms, quadratic forms, covariances, and coordinates before operations
that otherwise overflow or underflow. Positive covariance rescaling is allowed
only where the statistic/density is invariant to it. Do not independently
standardize coordinate axes or regressor columns without preserving the fitted
linear map and coefficient units.

Use stable differences, centered accumulation, expm1/log1p where appropriate,
log determinants, triangular solves, and log-sum-exp density accumulation.
Do not form an explicit covariance inverse merely to compute a quadratic form.
Solving once for a reusable factor is permitted. Projector notation does not
require an n by n projector allocation.

For a retained SPUR contrast with R'1=0, small positive decay must not be
evaluated by subtracting two nearly constant projected matrices. A stable
representative of the projected exponential covariance is

$$
R^\top\left[\frac{\operatorname{expm1}(-c d_{ij})}{2c}\right]_{ij}R
=-\tfrac12 R^\top[\,d_{ij}\phi(c d_{ij})\,]_{ij}R,
\qquad \phi(x)=-\operatorname{expm1}(-x)/x,\quad\phi(0)=1.
$$

It is proportional to R'K(c)R and has the Brownian limit at c=0. Use this
or equivalent stable arithmetic before trace normalization. The same
constant-annihilating `expm1(-c d)` form MAY be used for the centered SCPC
basis matrix $H_0=C_mK(c_0)C_m$ (and local centered landmark matrices) when
$\mathrm{avc}$ is close to one so small decay rates do not lose precision to
constant cancellation. The subtraction of the constant kernel is valid only
for a contrast annihilating constants; do not apply it to an unprojected SCPC
numerator. At short correlation scales, the unscaled kernel can be the more
stable representative. Any branch must preserve the same covariance shape and
pass the stated accuracy checks.

Let epsilon=2^-52. Symmetry and PSD checks MUST be scale-aware. For a p by p
projected covariance, the default entrywise symmetry allowance is
32 p epsilon times its largest absolute diagonal scale. A zero scale requires
an exact zero matrix or a separate documented construction-error bound.
Symmetrization within an assessed roundoff allowance is permitted; symmetrizing
a materially asymmetric input is not a repair.

PSD validation must check off-diagonal consistency, including a zero-variance
coordinate's cross-covariances. An eigenvalue or Schur-complement tolerance
must be normalized to the relevant block scale; one global scale MUST NOT hide
a materially indefinite small-variance block. Tiny negative eigenvalues
consistent with a documented backward-error estimate may be set to zero and
that adjustment must be available in diagnostics. Larger negativity is an
error. **[D]**

SCPC admits valid singular Gaussian covariance laws when their ratio statistic
is defined; see scpc.md for atom/tie behavior. SPUR angular densities and
Cholesky factors require positive definiteness on the retained contrast space.
Do not add jitter, a ridge, or silently remove requested contrasts to obtain
a factorization. Zero/failed response projections and unsupported singular
angular laws require a clear error.

### 2.1 Rank and eigenvalue blocks

The following dimension decisions are part of this profile:

| Quantity | Rule |
| --- | --- |
| Fitted-design estimability | Scale columns safely, use a rank-revealing decomposition, and compare singular values relative to the leading singular value at 1e-10 |
| Leading spatial eigenvalue eligibility | Require lambda_max>1.6e-11 for SCPC and the normalized Nyström pooled operator; require lambda_max>8 epsilon n for the unit-diameter SPUR Brownian construction |
| Spatial positive eigenvalues | Retain lambda > 1e-12 lambda_max, provided the leading value is above the applicable construction floor |
| SCPC eigenvalue-block gap | A block ends only at a gap greater than max(1.6e-11, 1.6e-11 abs(lambda_max)) for the unscaled centered K matrix, including each local Nyström landmark covariance |
| Nyström pooled eigenvalue-block gap | The same max(1.6e-11, 1.6e-11 abs(lambda_max)) rule on H_pool=EE' as defined in scpc.md, after each sample has unit Frobenius norm; use its own leading eigenvalue |
| SPUR eigenvalue-block gap | A block ends only at a gap greater than max(8 epsilon n, 1.6e-11 abs(lambda_max)) for the unit-diameter Brownian construction |
| Brownian whitening rank | Strict lambda > 1e-10 at unit diameter; defined further in half-life-and-transforms.md |
| Covariance shape equivalence | Frobenius norm of the difference of trace-normalized matrices <=1e-12, after finite positive traces and valid factors are established |

These thresholds are numerical support decisions, not scientific significance
levels. A user-specified q MUST end at a complete positive block. Automatic
selection MUST expand a partial eigensolve far enough to determine the retained
block boundary or return an error. Do not break a numerical tie arbitrarily to
obtain a requested dimension.

The Frobenius shape comparison is invariant under a common orthogonal change
of contrast coordinates. An entrywise comparison is not a substitute. Shape
closeness is a profile decision, not a bound on tail-probability differences
for ill-conditioned laws.

A partial eigensolve must verify residuals and orthogonality and resolve the
boundary gap before being accepted. Using a dense solve is an allowed fallback
when its memory requirement passes preflight; silently returning an unchecked
partial subspace is not.

The constant-span check uses the actual full design/FE projector. With Q an
orthonormal representation of its rank-r column space, require

$$
\|(I-QQ^\top)\mathbf 1/\sqrt n\|_2
\le 64\epsilon\max(1,r).
$$

An equivalent implicit projector is permitted. Verify its numerical accuracy;
if it cannot distinguish that condition, return a projection-precision error
rather than assuming a constant span.

For each **SPUR** response, form y_c=y-mean(y) safely by first subtracting a reference
entry and accumulating the shifted mean accurately; rescale by a common
power of two if necessary. A zero centered norm is a degeneracy error. For
residual SPUR let u=M_X y_c and require
||u||_2 > 1e-12 ||y_c||_2. For levels use u=y_c. For the spatial contrast
z=R' u require ||z||_2 > 1e-12 ||u||_2. These are numerical support thresholds,
not scientific model tests. Do not compare a small contrast with a huge
uncentered response offset. If construction error prevents resolving these
inequalities, report precision failure rather than pretending to distinguish
zero from nonzero. **[D]**

SCPC uses its fitted influence instead of this response-centering eligibility
rule. Let h_c be the centered site influence h_j and z=R' h_c. A zero h_c
or ||z||_2<=1e-12 ||h_c||_2 is a degenerate variance projection. For the
conditional denominator let V_pre=S_j' Xi_j R and
V_post=S_j' M' Xi_j R as defined in scpc.md. If V_pre has zero Frobenius norm,
or ||V_post||_F<=1e-12 ||V_pre||_F, report conditional-projection degeneracy.
This test is invariant under orthogonal rotation of the retained directions;
maxima of individual column norms are not substituted. Apply any common
rescaling consistently before comparing norms. A materially unresolved
construction error is a precision error. A constant outcome in a supported
no-intercept SCPC regression is not rejected merely for being constant.

The weight-norm check does not replace covariance-specific checks. If a required
SCPC covariance has exactly positive numerator variance and an identically
zero denominator block, follow scpc.md's calibration-degeneracy error before
starting a finite-cutoff solve. Co-location group cancellation can cause this
even when the weight-norm check passes. Distinguish algebraically established
zero from a numerically unresolved small variance.

### 2.2 Basis conventions and invariance

SCPC calibrations depend on the retained subspace, not its chosen orthonormal
orientation. SPUR statistics have the same distribution under orthogonal basis
rotation, but fixed simulation draws make orientation relevant to realized
results. **[I]**

Canonical coordinate ordering is used for geometry and landmark sampling.
Within a fixed implementation and numerical environment, use a deterministic
basis convention and reproduce seeded calls. No bespoke projector-anchor
construction or unique orientation of a tied eigenspace is required. Validate
the complete retained projector, eigen residuals and orthogonality; do not
use the tested response to orient or select directions. **[D]**

Identical supplied base Gaussian draws compare realized host outputs only when
their retained bases and coloring conventions agree. For a deliberate basis
change R_new=R_old O, an exact coupled comparison uses colored contrasts
G_new=O' G_old. With Cholesky factors for one covariance Omega, equivalent base
draws are Z_new=L_new^{-1} O' L_old Z_old, not generally O' Z_old. The map is
orthogonal for positive-definite covariances, but because Cholesky factors are
not rotation-equivariant, this map depends on Omega; multi-covariance
operations (such as the I(0) composite mesh or half-life search) therefore
require aligned retained bases R for identical finite-draw comparisons, or
distributional invariance checks otherwise. **[I]**

Repeated indistinguishable coordinate rows do not create an intrinsic order
among themselves. Do not invent persistent row identifiers, report an error
merely because an eigenbasis orientation is non-unique, or promise unchanged
fixed-seed ranks under every row permutation/backend change.

### 2.3 SCPC candidate expansion and landmark sampler

The initial working dimension is 10 for avc>=0.05, 20 for
0.01<=avc<0.05, 60 for 0.005<=avc<0.01, and 120 otherwise in the supported
domain. The dimension limit is m-1 in exact mode and min(m,landmarks)-1 in
Nyström mode. Cap each requested prefix at that limit.

If a prefix splits a positive block, or its largest candidate wins while more
positive directions remain available, increase the prefix to
min(limit,k+ceil(k/2)) and recompute candidates. A single landmark sample
does not disable this rule. Stop at an interior winner, exhaustion of positive
directions or the declared dimension limit, provided the selected prefix ends
at a complete block. If a required block cannot be resolved within the limit,
report the limitation instead of splitting it. This bounded search is an
efficiency heuristic. Geometry, decay limits, meshes and landmark subsets are
reused throughout expansion.

Sample uniform landmark subsets without replacement using a maintained,
versioned RNG implementation and its unbiased bounded-integer/subset sampler.
Do not prescribe the author replication's custom 32-bit LCG or retain a
legacy-bitstream inference mode. Record generator, sampler, version and seed;
use default seed 1. A fixed implementation must reproduce its subsets. **[D]**

Use the canonical site order and retain each selected subset throughout q
expansion. Subsets may overlap across samples; an unusable subset is reported,
not silently replaced by an outcome-dependent sample. Reconstructing the same
subsets must not consume a different stream. When ell=m use the complete-sample
shortcut in scpc.md. Cross-implementation approximation comparisons must align
the actual subsets, not merely assume equal seeds imply equal samples.

## 3. Common decay calibration

For m declared sites define

$$
\rho(c)=\frac{1}{m(m-1)}\sum_{i\ne j}\exp(-c d_{ij}).
$$

Solve rho(c)=rho_target only when f<rho_target<1 and there is a positive
distance. The SCPC default is avc=0.03 and its public domain is f<avc<1.
There is no additional inherited 0.001/0.99 restriction. Precision, rank and
resource checks determine whether a particular mathematically admissible
request can be computed. **[D]** An all-co-located geometry or unattainable
target is a support error, not a reason to change the requested correlation bound.

Evaluate the positive-pair contribution separately near f, and use stable
exponential differences near correlation one. A safeguarded monotone
bracket/Newton/bisection solver is the reference profile:

- At most 1,024 bracket expansions in either direction and 128 root iterations.
- Preferred correlation residual <=1e-12 times
  max(min(rho_target-f,1-rho_target),1e-15).
- A certified sign bracket of relative width <=1e-12 is also an acceptable
  stopping criterion; return which criterion was satisfied and the bracket.
- If neither criterion is satisfied, return a calibration accuracy error.
  A looser residual fallback MUST NOT be silently labelled converged.

The bracket is a bound for the numerically evaluated monotone function. It is
not an independent theorem about exact real-arithmetic evaluation. Detect
stagnation, nonfinite rates and failure to separate target from the floor.
Report the attained residual and bracket in diagnostic metadata.

Geometry summaries, f, positive minimum distance, and requested decay limits
MUST be computed once per geometry/parameter request and reused through q
expansion and coefficient calibration. Nyström uses all full-site pairs for
these quantities. Streaming pairs saves memory, not their quadratic time.

## 4. Shared covariance mesh

One bounded decay-rate mesh generator serves both SCPC (at caller $\rho_0=\mathrm{avc}$)
and I(0) (at $\rho_0=0.03$), solving only the two boundary rates rather than
running a separate root solve for every interior point. **[D]**

### 4.1 Interior and finite endpoint

For an attainable target $\rho_0\in(f,1)$, let $c_0$ solve $\rho(c_0)=\rho_0$
and let $c_{\max}$ solve

$$
\rho(c_{\max})=
f+\min\{(1-f)10^{-5},\ 0.1(\rho_0-f)\}.
$$

cmax is the upper endpoint of the **interior** mesh; it is not the largest
finite rate evaluated after tail augmentation. Report both cmax and c_tail.
The formula ensures a lower correlation target than $\rho_0$ when representable;
if a strictly positive span cannot be represented, return a clear calibration
error. **[D]**

The interior set is c0, all c0*1.05^k strictly between c0 and cmax, and cmax.
Use log-domain generation when necessary and include the actual endpoints
exactly once. Do not overshoot cmax merely to force a minimum point count.
**The interior mesh MUST NOT exceed 1,024 steps** (requiring
$\lfloor(\log c_{\max}-\log c_0)/\log 1.05\rfloor\le 1{,}024$); if an extreme
$\rho_0$ near $1$ or multiscale geometry exceeds this budget, fail at preflight
with an explicit mesh-budget error naming the span and limit rather than
coarsening silently. The 1.05 spacing is a retained numerical choice, not an
accuracy certificate.

### 4.2 Shared bounded tail construction

Let d_min be the smallest strictly positive full-site distance in the same
working units as c. From $a=c_{\max}>0$, define

$$
b=\max\{a,\ 36/d_{\min}\},\qquad
N=\left\lceil\frac{\log b-\log a}{\log 1.2}\right\rceil.
$$

If b=a, the tail contributes no new finite points. Otherwise include

$$
c_k=\exp\{\log a+(k/N)(\log b-\log a)\},\quad k=1,\ldots,N,
$$

with c_N assigned exactly b. Require finite, strictly increasing representable
rates. **N must not exceed 256.** If b is unrepresentable, N exceeds 256, or
progress is lost, return an explicit tail-calibration/resource error naming
the distance range and required budget. Do not cap b silently, coarsen the
grid, or switch to an undocumented dense-window fallback.

Add every tail point and K(infinity), using the true co-location block
covariance. The finite mesh depends only on geometry and $\rho_0$, and MUST NOT
depend on the tested response. For SCPC, use this identical set for:

1. Every q candidate's location-only objective.
2. Every reported location and conditional critical value.
3. Every coefficient p-value and reported confidence level.

The full-site projected family, including the same limit, also applies when
the basis is Nyström. A cached covariance may be reused or tiled rather than
stored, but it must not be omitted from any of those calculations.

At 36/d_min all positive-distance kernel entries are at most exp(-36).
That is an entrywise kernel observation. It does **not** bound projected
probability error under ill-conditioning or prove that a finite interior
maximum has been located. **[I, D]**

### 4.3 Reporting and refinement

Return calibration_status=finite_grid, the actual finite rate mesh, units,
c0/cmax/c_tail, f, and inclusion of the exact limit. Record the maximizing
evaluated covariance for each reported probability/cutoff. The standard
profile does not promise a continuous calibration error bound.

A diagnostic sensitivity run may add midpoint/tail points and report changes.
It must not silently change the production family, select points from the
response without declaration, or be presented as a global certificate.
Promoting a different production mesh requires a new profile record.

## 5. I(0) covariance mesh

Use the SPUR unit-diameter geometry. If f>0.03 the specified null family is
infeasible: return the error defined in spur.md. If f=0.03, the null consists
only of K(infinity); evaluate that exact limit and do not construct a finite
rate grid or tail.

For f<0.03, construct the I(0) nuisance mesh with the **same §4 bounded decay
mesh** at $\rho_0=0.03$: solve only $c_{0.03}$ and $c_{\max}$ (for
$\rho(c_{\max})=f+\min\{(1-f)10^{-5},\ 0.1(0.03-f)\}$), generate the 1.05
interior mesh (`<= 1,024` steps), append the 1.20 tail up to
$b=\max(c_{\max},36/d_{\min})$ (`<= 256` steps), and include K(infinity).
This unifies SCPC and I(0) on one two-solve decay-mesh implementation instead
of maintaining a separate ~79-solve hybrid correlation grid and deduplication
pass, while keeping interior rate ratios at most 1.05 and enforcing the same
preflight budget errors. **[D]**

The identical nuisance set MUST determine I(0) rank p-values and every
reported cutoff. Tuning has its own fixed reference-null experiment as
defined in spur.md; the reference tuning null must not be confused with this
composite inference family. Report the actual rates, achieved average
correlations, limit, maximizing points, and finite_grid status.

## 6. Gaussian probabilities and SCPC critical values

Use a tail-stable Gaussian quadratic-form calculation. The current angular
representation and its complementary form are the reference approach; a
validated mathematically equivalent solver is permitted. The solver must
distinguish strict from inclusive tails for genuine atoms as defined in scpc.md.
An unsupported atom or numerically unresolved singular ratio requires a
descriptive error rather than an invalid strict-tail p-value. Check the exact
zero-denominator/positive-numerator case before root finding: its cutoff is
mathematically infinite, and this profile returns the explicit
calibration-degeneracy error in scpc.md. It is not a failure to bracket an
otherwise finite root. A zero numerator and denominator law is undefined.

For ordinary probabilities, the default requested absolute integration accuracy
is 1e-11. For a tail or complement smaller than 1e-8, request relative accuracy
1e-8 for that smaller probability; if it lies outside the representable range,
report the numerical limit instead of inventing a positive floor.
These are numerical error criteria, not Monte Carlo standard errors.
The relative criterion applies to the directly evaluated small tail or
complement. A binary64 public p-value may round to one when its accurately
computed complement is too small to survive subtraction. Do not claim that
forming 1-p recovers that complement at the same relative accuracy; retain
its directly computed value in diagnostics when needed to explain the limit.

The reference adaptive integrator has maximum recursion depth 20. At the cap
it MUST either have satisfied its declared error criterion or return an
accuracy failure. Record an estimated integration error where the estimator
provides one, and distinguish an estimate from a rigorous bound. Clamping a
probability to [0,1] is permitted only for an excursion within its assessed
roundoff/error allowance; larger excursions are errors.

Critical-value tolerance is 0.001 in ordinary t units by default. A caller may
request a smaller positive tolerance within supported precision. Bracket the
smallest cutoff satisfying the required evaluated-family tail criterion:

- A Student-t value may be an initial guess, never a lower bound.
- Search downward as well as upward when needed.
- Preserve a rejecting lower and non-rejecting upper bracket, including
  separately handled atom boundaries.
- Return the upper endpoint only after its width is within the requested
  tolerance and the final upper criterion is checked over the complete family.
- Record the final bracket and maximum achieved tail. Iteration/precision
  stagnation without a sufficiently narrow valid bracket is an accuracy error.

For continuous laws, p/interval disagreement due only to upper-bracket rounding
must be confined to that returned t bracket, or to a case where the assessed
probability-integration uncertainty overlaps alpha. Use error-aware bracket
comparisons (upper estimated tail plus its assessed error <=alpha; lower
estimated tail minus its assessed error >alpha), refining the integral when
needed to distinguish a bracket endpoint. Do not label an error estimate a
rigorous bound. A larger discrepancy outside both uncertainty allowances, or
a different nuisance family for p and cutoff, is a conformance failure.
Exact equality uses the closed-interval/strict-cutoff convention in scpc.md.

## 7. Simulations, tuning, and half-life quadrature

The default SPUR inference block has B=100000 independent standard Gaussian
q-vectors. Diagnostics use an independent tuning block with its own count
B_tune, also defaulting to 100000; the counts need not match. Require B>=1
and B_tune>=19, so the 5% tuning test can have a finite rank cutoff.
Half-life uses one inference block because its alternative is geometry-only.
Reuse the same base Gaussian vectors across parameter points and recolor them
separately by each null covariance's specified Cholesky factor. A radial rescaling
of each column does not change the homogeneous statistics. **[D, I]**

Accept user-supplied blocks only after checking dimensions, finite entries, and
strictly positive norms. State the required probabilistic generation and
independence contract; those properties cannot be established by a norm check.
User-supplied arbitrary deterministic quadrature points do not acquire
finite-simulation rank validity merely by passing input validation.

The point-null rank guarantee concerns the joint sampling law of the observed
statistic and its independent inference block, conditional on the eligible
geometry/design and independent tuning result where applicable. A fixed
realized draw matrix supports reproducibility; it is not a conditional size
certificate for every possible realized simulation matrix.

Minimum attainable p is 1/(B+1). This is rank resolution, not probability
accuracy: a fixed point-null tail near 0.05 has binomial Monte Carlo standard
error about 0.00069 at B=100000. No such scalar error bar certifies the maximum
over a composite family or the completeness of a confidence-set search.
Rank/tie and cutoff definitions are in spur.md
and half-life-and-transforms.md. Finite-simulation inference is a numerical
profile in its own right; do not label its p-values exact deterministic Gaussian
probabilities. Common seeds across different host RNGs do not imply common
draws.

The ten-node geometry-dependent composite half-life quadrature and complete
bounded search recipe are defined in half-life-and-transforms.md, governed by
the single 64-panel (640-node) short-scale preflight cap (which also bounds
generated short-scale search probes to at most 441). That chapter owns their
constants. This profile retains them because independent convergence
checks supported ten-node integration; it does not increase orders or rounds
without evidence of a material remaining error.

The bounded reference-tuning search is fully defined in spur.md. It uses
unit-trace covariance shapes, a strength parameter in [0,1], exact endpoint
laws, and the same rank-cutoff helper as inference. There are no type-7
quantiles, unbounded multiplier searches, author-code starting constants or
floating-point-underflow stopping rules. Failure to reach the target is
reported with achieved estimated power and status; it must not be advertised
as final-test power or global optimum.

## 8. Efficiency and resource contract

Do not repeat geometry-only decay calibration during q expansion. Reuse
orthogonal projections, covariance factors, the half-life mixture, observed
mixture evaluations, and repeated point evaluations within one call, with
dependencies including sample, design, geometry, q, numerical profile, and
draw blocks as applicable. A persistent global cache is not required.

Tile projected covariances and simulation work. Before allocation, check
dimension products, platform integer limits, anticipated workspace and caller/
environment resource limits. Report an actionable resource error instead of
overflow, process abort, or silently reducing q/B/landmarks/search resolution.
An automatically chosen Nyström basis is declared as such and still uses the
full-site covariance mesh; it is not a secret fallback after failed exact
inference.

The SCPC finite meshes are constructed before outcome inference and are finite
by the above rules. Half-life's fixed probes, refinement rounds and boundary
budget are its explicit algorithmic limits. A host may enforce a tighter
resource budget, but must reject before work or return an incomplete diagnostic
result; it must not report normal completed inference.

Performance changes should preserve this contract. In particular:

- Full-pair streaming is O(m^2) time even if pair storage is O(m).
- Nyström's auto threshold is a memory policy, not a universal speed claim.
- More uniform grid points do not by themselves justify a stronger guarantee.
- Solver and search benchmarks must record n/m, q, B, profile, seeds, threads,
  cold/warm runs, memory and outputs.

## 9. Optional deterministic SPUR solver

For a fixed covariance and ratio of quadratic forms, a general validated
Davies/Imhof calculation is an allowed future numerical profile. It is not a
required 0.6.0 production rewrite or an unmeasured speed claim. Promotion
requires general indefinite-form accuracy/failure handling, independent
probability checks, realistic timing/memory evidence and clear p-value/cutoff
semantics. The specialized SCPC solver is not automatically sufficient.
See [Davies](sources.md#davies).

Changing from finite rank inference to deterministic probabilities requires
explicit profile/version documentation. It does not solve continuous nuisance
optimization or the half-life mixture-ratio distribution. No legacy simulation
path is required solely to reproduce an obsolete approximation once an
intentional method/profile change has been adopted and verified.
