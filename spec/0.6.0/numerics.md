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
is defined; see scpc.md for atom/tie behavior. While a unique symmetric
positive-semidefinite square root $\Omega^{1/2}$ exists mathematically for
singular PSD matrices as well, the adopted SPUR and half-life angular-density
and inverse-covariance quadratic-form experiment requires positive definiteness
on the retained contrast space. Do not add jitter, a ridge, or silently remove
requested contrasts to obtain a factorization. Zero/failed response projections
and unsupported singular angular laws require a clear error.

### 2.1 Rank and eigenvalue blocks

The following dimension decisions are part of this profile. For any
$m_{\mathrm{op}}\times m_{\mathrm{op}}$ spatial operator (unscaled centered SCPC
$H_0$ with $m_{\mathrm{op}}=m$, local Nyström $H_r$ with $m_{\mathrm{op}}=\ell$,
Frobenius-normalized pooled Nyström $H_{\mathrm{pool}}$ with $m_{\mathrm{op}}=m$,
or unit-diameter SPUR Brownian $B_X$ with $m_{\mathrm{op}}=n$), define the
unified construction floor

$$
\tau_{\mathrm{floor}}=64\epsilon m_{\mathrm{op}},
\qquad \epsilon=2^{-52}.
$$

| Quantity | Rule |
| --- | --- |
| Fitted-design estimability | Scale columns safely, use a rank-revealing decomposition, and compare singular values relative to the leading singular value at 1e-10 |
| Spatial positive eigenvalues | Retain $\lambda>\max(\tau_{\mathrm{floor}},\ 10^{-12}\lambda_{\max})$, using the operator's own size $m_{\mathrm{op}}$ and leading value $\lambda_{\max}$ (reporting $\lambda_{\max}$ and $\tau_{\mathrm{floor}}$ if no eligible complete block exists) |
| Spatial eigenvalue-block gap | A block ends only at a gap $\lambda_k-\lambda_{k+1}>\max(\tau_{\mathrm{floor}},\ 10^{-11}\lambda_{\max})$, using each operator's own $m_{\mathrm{op}}$ and $\lambda_{\max}$ |
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

Identical supplied base Gaussian draws compare realized host outputs when their
retained bases agree. For any orthogonal change of retained basis
$R_{\mathrm{new}}=R_{\mathrm{old}}O$, SPUR and half-life color draws with the
unique symmetric positive-definite square root
$S(\Omega)=\Omega^{1/2}=V\Lambda^{1/2}V^\top$. Because
$(O^\top\Omega O)^{1/2}=O^\top\Omega^{1/2}O$, rotating the base draws once by
$E_{\mathrm{new}}=O^\top E_{\mathrm{old}}$ gives an exact coupled comparison
$S(\Omega_{\mathrm{new}})E_{\mathrm{new}}=O^\top S(\Omega_{\mathrm{old}})E_{\mathrm{old}}$
simultaneously for **every** covariance $\Omega$. Consequently, both single-covariance
and multi-covariance operations (I(1), the I(0) composite maximum, multi-candidate
tuning, and the jointly searched half-life confidence set) have exact joint
finite-$B$ distributional invariance under orthogonal basis rotations without
any $\Omega$-dependent draw transformation. Implementations MAY reuse the same
symmetric eigendecomposition $\Omega=V\Lambda V^\top$ already needed for
$\Omega^{1/2}$ to evaluate quadratic forms $\|\Lambda^{-1/2}V^\top z\|^2$ and
$\log|\Omega|=\sum_i\log\lambda_i$, or use Cholesky/triangular solves,
provided accuracy checks on ill-conditioned covariances pass. **[I, D]**

Repeated indistinguishable coordinate rows do not create an intrinsic order
among themselves. Do not invent persistent row identifiers, report an error
merely because an eigenbasis orientation is non-unique, or promise unchanged
fixed-seed ranks under every row permutation/backend change.

### 2.3 SCPC candidate expansion and landmark sampler

The initial working dimension is $\min(\mathrm{limit},20)$ across all supported
$\mathrm{avc}$ values, where $\mathrm{limit}$ is $m-1$ in exact mode and
$\min(m,\mathrm{landmarks})-1$ in Nyström mode.

If a prefix splits a positive block, or its largest candidate wins while more
positive directions remain available, double the prefix to $\min(\mathrm{limit},2k)$
and recompute candidates. A single landmark sample
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
expansion. Subsets may overlap across samples; any unusable subset (zero
positive rank or zero/nonfinite embedded norm) MUST fail clearly with a
recommendation to increase landmarks or use exact geometry, rather than
skipping subsets or resampling subsets automatically. Reconstructing the same
subsets must not consume a different stream. When ell=m use the complete-sample
shortcut in scpc.md. Cross-implementation approximation comparisons align
subsets by sharing the same compiled core RNG and seed, or by explicitly
supplying/comparing the actual landmark indices.

## 3. Common decay calibration

Covariance methods normalize pairwise distances internally by the positive
sample diameter $\Delta$ ($d_{ij}=D_{ij}/\Delta\in[0,1]$, $d_{\max}=1$) and
convert reported decay rates at the API boundary via $c_{\mathrm{physical}}=c/\Delta$.
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

Because $0<d_{\min}\le d_{ij}\le d_{\max}\le1$ on positive-distance pairs,
defining

$$
a=\frac{\rho_{\mathrm{target}}-f}{1-f},\qquad
L=\begin{cases}
-\log(a), & a\le\tfrac12,\\
-\operatorname{log1p}\!\left(\dfrac{\rho_{\mathrm{target}}-1}{1-f}\right), & a>\tfrac12,
\end{cases}
$$

gives the exact mathematical analytic bounds $L/d_{\max}\le c_0\le L/d_{\min}$
(which under unit-diameter normalization $d_{\max}=1$ simplify to
$L\le c_0\le L/d_{\min}$). Using $-\log(a)$ for $a\le\tfrac12$ avoids
cancellation near the co-location floor $f$, while $-\operatorname{log1p}$ for
$a>\tfrac12$ avoids cancellation near correlation $1$.

When evaluating residuals $\rho(c)-\rho_{\mathrm{target}}$:

- Near the co-location floor $f$, evaluate the positive-distance contribution
  separately as $\frac{1}{m(m-1)}\sum_{d_{ij}>0}\exp(-c d_{ij})-(\rho_{\mathrm{target}}-f)$
  rather than adding $f$ and subtracting $\rho_{\mathrm{target}}$.
- Near correlation $1$, evaluate correlation deficits using
  $\operatorname{expm1}(-c d_{ij})$ rather than forming $1-\exp(-c d_{ij})$.

If $d_{\min}=d_{\max}$, $c_0=L/d_{\max}$ holds in closed form and MUST still be
checked against the residual accuracy criterion below. Otherwise validate that
both analytic bracket endpoints are finite and representable, check the residual
criterion at both endpoints (in case floating-point roundoff places the root on
an endpoint), and run one safeguarded Newton/Brent/bisection scalar solver
inside that bracket using the following reference profile:

- No exploratory upward/downward bracket-expansion loop is required; at most 128
  scalar solver iterations inside $[L/d_{\max},L/d_{\min}]$.
- Preferred correlation residual <=1e-12 times
  max(min(rho_target-f,1-rho_target),1e-15).
- A certified sign bracket of relative width <=1e-12 is also an acceptable
  stopping criterion; return which criterion was satisfied and the bracket.
- If neither criterion is satisfied or bracket endpoints cannot be represented,
  return a calibration accuracy error. A looser residual fallback MUST NOT be
  silently labelled converged.

The analytic bounds $[L/d_{\max},L/d_{\min}]$ hold identically in exact real
arithmetic; by contrast, the rounded, numerically verified sign bracket returned
by the solver bounds only the numerically evaluated monotone function and is not
an independent theorem about exact real-arithmetic evaluation. Detect
stagnation, nonfinite rates and failure to separate target from the floor.
Report the attained residual and final numerical bracket in diagnostic metadata.

Geometry summaries, f, positive minimum distance, and requested decay limits
MUST be computed once per geometry/parameter request and reused through q
expansion and coefficient calibration. Nyström uses all full-site pairs for
these quantities. Streaming pairs saves memory, not their quadratic time.

## 4. Shared covariance mesh

One single-solve uniform log-rate mesh generator serves both SCPC (at caller
$\rho_0=\mathrm{avc}$) and I(0) (at $\rho_0=0.03$), eliminating separate interior
vs. tail zones and second root solves. **[D]**

### 4.1 Single bounded log-rate construction

For an attainable target $\rho_0\in(f,1)$, solve $\rho(c_0)=\rho_0$ for $c_0>0$.
Let $d_{\min}$ be the smallest strictly positive full-site distance in the same
working units as $c$, and define

$$
c_{\mathrm{end}}=\max\{c_0,\ 36/d_{\min}\}.
$$

If $c_{\mathrm{end}}=c_0$, the finite mesh is $\{c_0\}$ ($N=0$, 1 finite point).
Otherwise ($c_{\mathrm{end}}>c_0$), set

$$
N=\left\lceil\frac{\log c_{\mathrm{end}}-\log c_0}{\log 1.10}\right\rceil,\qquad
c_k=\exp\!\left(\log c_0+\frac{k}{N}(\log c_{\mathrm{end}}-\log c_0)\right),\quad k=0,1,\ldots,N,
$$

with exact endpoints $c_0$ and $c_N=c_{\mathrm{end}}$ ($N+1$ finite points).
Require finite, strictly increasing representable rates. **$N$ MUST NOT exceed
512** (covering a $1.10^{512}\approx1.5\times10^{21}$ decay span). If
$c_{\mathrm{end}}$ is unrepresentable, $N>512$, or numerical progress is lost,
reject at preflight with an actionable mesh-size error reporting $c_0$,
$c_{\mathrm{end}}$, $d_{\min}$, required $N$, and cap `512`. Do not cap
$c_{\mathrm{end}}$ silently, coarsen the grid, or switch to an undocumented
dense-window fallback.

Add the exact co-location limit $K(\infty)$ ($N+2$ total covariance points when
$c_{\mathrm{end}}>c_0$), using the true co-location block covariance. The mesh
depends only on geometry and $\rho_0$, and MUST NOT depend on the tested
response. For SCPC, use this identical set for:

1. Every q candidate's location-only objective.
2. Every reported location and conditional critical value.
3. Every coefficient p-value and reported confidence level.

The full-site projected family, including the same limit, also applies when
the basis is Nyström. A cached covariance may be reused or tiled rather than
stored, but it must not be omitted from any of those calculations.

At $36/d_{\min}$ all positive-distance kernel entries are at most $\exp(-36)\approx2.3\times10^{-16}$.
That is an entrywise kernel observation. It does **not** bound projected
probability error under ill-conditioning or prove that a finite interior
maximum has been located. The `1.10` step ratio is a retained practical
numerical choice, not a continuous accuracy certificate. **[I, D]**

### 4.2 Reporting and refinement

Return `calibration_status=finite_grid`, the actual finite rate mesh, units,
$c_0$, $c_{\mathrm{end}}$, $f$, and inclusion of the exact limit $K(\infty)$.
Record the maximizing evaluated covariance for each reported probability/cutoff.
The standard profile does not promise a continuous calibration error bound.

A diagnostic sensitivity run may evaluate a finer grid (e.g., ratio `1.05` or
closer) and report changes. It must not silently change the production family,
select points from the response without declaration, or be presented as a global
certificate. Promoting a different production mesh requires a new profile record.

## 5. I(0) covariance mesh

Use the SPUR unit-diameter geometry. If $f>0.03$ the specified null family is
infeasible: return the error defined in spur.md. If $f=0.03$, the null consists
only of $K(\infty)$; evaluate that exact limit and do not construct a finite
rate grid.

For $f<0.03$, construct the I(0) nuisance mesh with the **same §4 single
log-rate mesh** at $\rho_0=0.03$: solve $c_{0.03}$, set
$c_{\mathrm{end}}=\max(c_{0.03},36/d_{\min})$, generate the single
log-rate grid at step ratio at most `1.10` (`N <= 512`), and include
$K(\infty)$. This unifies SCPC and I(0) on one single-solve decay-mesh
generator instead of maintaining a separate ~79-solve hybrid grid or a stitched
two-zone mesh, while enforcing the same preflight budget error and full I(0)
validation checks in verification.md. **[D]**

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

- Initialize from a simple positive starting value or previously computed
  cutoff, never imposing a Student-$t$ lower floor (no Student-$t$ quantile
  routine is required).
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
separately by each null covariance's unique symmetric positive-definite square
root $\Omega^{1/2}$ (so orthogonal basis rotations preserve the joint finite-$B$
simulation distribution across all parameter points). A radial rescaling of
each column does not change the homogeneous statistics. **[D, I]**

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
probabilities. Ordinary seeded draws are generated inside the shared core via one
maintained RNG, so hosts sharing the same compiled core binary and seed produce
identical generated blocks; cross-version or custom-generator comparisons use
explicitly supplied draws.

The ten-node geometry-dependent composite half-life quadrature and bounded
uniform log-h search recipe (including their panel cap, initial grid, and
bisection rules) are defined normatively in half-life-and-transforms.md, which
serves as their single normative home. This profile retains ten-node quadrature
because independent convergence checks supported it; it does not increase
integration orders without evidence of a material remaining error.

The bounded reference-tuning search is fully defined in spur.md. It uses
unit-trace covariance shapes, a strength parameter in [0,1], exact endpoint
laws, symmetric square-root coloring, and the same rank-cutoff helper as
inference. There are no type-7 quantiles, unbounded multiplier searches,
author-code starting constants or floating-point-underflow stopping rules. Failure
to reach the target is reported with achieved estimated power and status; it must
not be advertised as final-test power or global optimum.

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
by the above rules. Half-life's uniform log-h grid and 20-step bisection are its
explicit algorithmic limits. A host may enforce a tighter resource budget and
fail with an actionable error at preflight or during execution. It MAY return an
explicitly incomplete diagnostic result; it must not report normal completed
inference.

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
