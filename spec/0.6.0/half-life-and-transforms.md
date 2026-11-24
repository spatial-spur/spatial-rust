# Half-life inference, spatial transformations, and the combined workflow

## 1. Status, provenance, and scope

**This chapter is a normative target for methods version 0.6.0. Publishing the
chapter does not establish that an existing binary conforms to it.** An
implementation may advertise conformance only after satisfying the applicable
[verification requirements](verification.md). In particular, the complete
evaluation trace, incomplete-search reporting, and conservative uninformative
result defined below are requirements, not claims about version 0.5.0. [D]

Provenance tags have the meanings defined in the [specification overview](README.md):
**[P]** published mathematics; **[R]** published tuning or applied recommendation;
**[A]** author-software convention; **[D]** package decision; **[I]** algebraic
consequence. Source identities and replication-file locations are recorded in
[sources](sources.md). Statements that combine a mathematical target with a
numerical implementation carry more than one tag.

Half-life inference concerns a scalar response under the canonical spatial
Gaussian covariance family. Its broader interpretation requires the spatial
asymptotic assumptions stated in the [SPUR chapter](spur.md); numerical rank,
finite input values, or successful matrix factorization do not establish those
assumptions. The transformation definitions below are separate from guarantees
for a subsequent regression. [P, D]

## 2. Half-life geometry, experiment, and eligibility

Let the analysis sample contain $n$ observations at locations $s_i$, with
distance matrix $D$, positive diameter
$d_{\max}=\max_{i,j}D_{ij}$, and
$M=I-\mathbf1\mathbf1^\top/n$. Distances, exact co-location, geographic arc
units, and input ordering follow [the common specification](README.md) and
[SPUR geometry](spur.md). Define normalized distances $d=D/d_{\max}$ and

$$
 B=-\tfrac12 M d M.
$$

Let $R_q$ have unit-norm orthogonal columns spanning the leading eligible
positive eigenspaces of $B$, and let $z=R_q^\top M y$. Eligibility, complete
eigenvalue blocks, basis conventions, and numerical rank follow [SPUR](spur.md)
and [numerics](numerics.md). The exact centered Brownian covariance implements
the published target; the author replication's high-correlation exponential proxy
is not part of this definition. [P, I, D]

The supported half-life operation requires an integer $2\le q<n$, sufficient
positive retained rank, a permitted eigenspace cutoff, and $q\le100$. The
default is $q=15$. The lower bound and ceiling define the supported numerical
operation; neither is a universal theoretical validity boundary. The ceiling
is a workload and validation policy, not a uniform quadrature error bound.
Retaining high frequencies or all available positive directions requires the
advisories specified in [SPUR](spur.md). [R, D]

Write normalized half-life as $h$, with physical half-life
$H=d_{\max}h$. For $0<h<\infty$, set $c=\log(2)/h$ and use the
projected covariance shape

$$
 \Omega(h)\ \propto\ R_q^\top K(c)R_q,
 \qquad K(c)_{ij}=\exp(-c d_{ij}).
$$

One paper-normalized representative is
$R_q^\top K(c)R_q/(2c)$. Positive scalar multiplication of any covariance
does not change the angular experiment. [P, I]

The endpoint shapes are

$$
 \Omega(0)=R_q^\top J_0R_q,\qquad
 (J_0)_{ij}=\mathbf1[d_{ij}=0],
 \qquad
 \Omega(\infty)=-\tfrac12 R_q^\top dR_q.
$$

Thus co-located observations retain their exact covariance block at zero
half-life. They are not silently replaced by independent observations or
jittered coordinates. With distinct locations, $\Omega(0)=I_q$. The
infinite-half-life experiment is evaluated as a limit, not by a finite
sentinel subsequently displayed as infinity. [I, D]

All required retained covariance shapes must admit the positive-definite
factorizations defined by [numerics](numerics.md). Invalid covariance geometry,
insufficient rank, a non-finite response, a constant response, or a numerically
zero projected vector $z$ produces an actionable error. An undefined
likelihood is not an empty confidence set. [D]

The public half-life operation is a levels-response operation. An extension to
residual half-life requires an explicit nuisance-design contract and is not
implied merely by support for residual I(0)/I(1) diagnostics. [D]

## 3. Angular likelihood and point-null tests

### 3.1 Mathematical target

For $z\ne0$, set $u=z/\lVert z\rVert$. Up to a common constant depending
only on $q$, the angular Gaussian density is

$$
 f_\Omega(u)=|\Omega|^{-1/2}
             (u^\top\Omega^{-1}u)^{-q/2}.
$$

The target alternative averages these **densities**, with a uniform measure
on normalized half-life:

$$
 m(u)=\int_0^1 f_{\Omega(h)}(u)\,dh,
 \qquad T_{h_0}(u)=\frac{m(u)}{f_{\Omega(h_0)}(u)}.
$$

It does not average covariance matrices and is not uniform in decay $c$.
Null values $h_0$ range over $[0,\infty]$; the alternative's support does
not truncate the null domain. In physical units the same mixture has measure
$dH/d_{\max}$ on $[0,d_{\max}]$. This is the equation-(20) target of
Müller and Watson (2024), §4.4, with constants that do not affect a calibrated
point-null test included in the ratio. [P, I]

The same homogeneous expression can be evaluated using unnormalized $z$
in numerator and denominator, since its common response-scale factor cancels.
Covariance scaling, omission of the uniform normalizer in physical units, or
omission of a null-specific positive determinant factor from both observed
and simulated scores must not change the point-null ranking. [I]

### 3.2 Finite simulation and acceptance

The public API must provide a point-null operation at any requested
$h_0\in[0,\infty]$, using the same geometry, quadrature profile, and supplied
or generated inference draws as the confidence-set operation. The operation
preserves a supplied finite $h_0$ and evaluates its finite-half-life kernel
directly; it must not first round $h_0/(1+h_0)$ to one and silently replace
that request by the infinite-half-life experiment. It
requires a finite confidence level $0<\mathrm{level}<1$, default 0.95,
and sets the nominal rejection size to $\alpha=1-\mathrm{level}$.
An adapter exposing percentages converts once to this normalized binary64
level, for example 95 to 0.95, and reports the normalized value used.
It returns the tested half-life and units, p-value, acceptance at the requested
level, and numerical status. A simulated point test additionally returns its
actual exceedance and draw counts. Analytic or conservative-policy results
identify their source and do not fabricate simulation counts. [D]

For a fixed null, generate $Z_b=L_{h_0}\xi_b$, where
$L_{h_0}$ is the positive-diagonal lower Cholesky factor of the chosen
positive covariance representative in the retained basis used by the operation, following
the fixed-draw convention in [SPUR](spur.md#62-ranks-cutoffs-and-equality).
The columns
$\xi_b$ obey the distribution and independence contract in [SPUR](spur.md).
Use the same deterministic statistic for the observation and each simulation.
With $B_{\mathrm{MC}}$ inference draws, define

$$
 N(h_0)=\#\{b:T_{h_0}(Z_b)\ge T_{h_0}(z)\},\qquad
 p_B(h_0)=\frac{1+N(h_0)}{B_{\mathrm{MC}}+1}.
$$

**Accept if $p_B(h_0)>\alpha$; reject if $p_B(h_0)\le\alpha$.**
For a binary64 requested level, implement the profile's attainable-boundary
convention by forming the integer $K=B_{\mathrm{MC}}-N(h_0)$ and comparing
the binary64 quotient $K/(B_{\mathrm{MC}}+1)$ directly with that level,
using ordinary correctly rounded binary64 count conversion and division:
accept exactly when it is strictly smaller. A level with the
same binary64 representation as that complementary rank boundary rejects;
do not form $1-\mathrm{level}$ to make the acceptance decision or add an
epsilon band. This is a numerical nominal-level convention, not a promise of
exact real-rational comparison for every arbitrarily supplied real level.
Integer arithmetic must be checked for overflow. Ties count as at least as
extreme. The default draw count is 100,000; draw identity and RNG metadata are
reported under the common simulation contract. [D]

Common inference draws are reused across null half-lives. This does not require
independence between tests at different $h_0$. The point-null rank argument
requires appropriate simulated directions independent of the observation;
finite positive column norms alone do not establish that statistical property.
Column normalization is allowed because the score is homogeneous. [I, D]

A fixed positive quadrature mixture defines a legitimate test statistic even
when it only approximates $m$. Its quadrature error concerns fidelity to the
specified mixture and its power properties; it is distinct from invalid null
simulation or an incomplete search for accepted parameters. No quadrature
adaptation may use the observed outcome to select a different statistic from
the one applied to simulated observations. [I, D]

### 3.3 Uninformative experiments

If the valid, nonzero projected experiment has covariance shapes that are
structurally proportional for all half-lives, no angular statistic identifies
half-life. Return $p(h_0)=1$ for every point-null request and the exact
full-domain set $[0,\infty]$, with an explicit structural-uninformative
status. Unsupported dimensions, such as $q=1$ in this profile, need not be
accepted to implement this behavior. [I, D]

The profile also permits a **conservative uninformative fallback** when its
geometry-only shape-equivalence criterion in [numerics](numerics.md) is met
over the required initial and mixture covariances. Return $p(h_0)=1$ and the
full domain, labelled `uninformative_within_profile`. This is a deliberate
conservative test/reporting policy. Equality on a finite mesh or within a
tolerance is not a proof that the continuum has identical covariance shapes.
Record the criterion and evaluated covariances that triggered the fallback.
Do not substitute it for a failed factorization, zero projected response,
interrupted computation, or exhausted numerical accuracy check. [D]

If the requested level exceeds
$B_{\mathrm{MC}}/(B_{\mathrm{MC}}+1)$, rejection is unattainable under the
rank rule and a valid experiment necessarily returns the full domain. Report
that reason separately from lack of spatial information. Equality does permit
rejection. Apply the same represented-boundary convention from §3.2 to this
comparison. The implementation may short-circuit this case after eligibility
validation instead of evaluating unnecessary simulation points. [I, D]

These full-domain shortcuts do not require simulation solely to populate
result fields. Record only work actually performed, with empty point-evaluation
traces and absent counts permitted when no point test was run. A covariance
checked for eligibility or shape equivalence is not a simulated point-null
evaluation. Record the reason when exact zero or infinity point tests are
bypassed. The unattainable-rank-level shortcut establishes acceptance over the
domain, not that every point-null p-value equals one; an independently requested
point test still returns its actual rank p-value unless an uninformative policy
applies. [I, D]

## 4. Retained mixture quadrature profile

Let $d_{\min}$ be the smallest strictly positive normalized inter-site
distance and let $\nu$ be the smallest positive normal binary64 number.
Set $a=\max(d_{\min}/32,\nu)$. Starting from $r=2^{-9}$, halve $r$
until $r\le a$. The panels are $[0,r],[r,2r],[2r,4r],\ldots$, ending
at 1. This profile supports **at most 64 panels**, or 640 mixture nodes.
Determine the required panel count from geometry before constructing factors
or simulations. Exceeding the cap is a resource/geometry error reporting
$d_{\min}$, the sample diameter, the required count and the limit; do not
merge near locations, clip the scale or substitute the zero-half-life law.
Each panel uses the standard **ten-point Gauss–Legendre** nodes and positive
integration weights transformed to that panel. The weights sum to one to
numerical accuracy. [D]

This preserves the current geometry-aware ten-node policy. The small first
panel resolves scales below the closest positive separation: for
$h\le d_{\min}/32$, every positive-distance kernel correlation is at most
$2^{-32}$. This observation motivates the panel boundary; it is not a
uniform bound on angular-mixture error for every projected experiment. [I, D]

Construct the nodes and weights from geometry only, once per prepared
operation. Cache their projected covariance factors and determinant terms.
Use stable factor solves and log accumulation, or an equivalent calculation
meeting [numerics](numerics.md). A finite positive half-life does **not**
require a representable intermediate decay $c=\log(2)/h$. For short scales,
evaluate $\exp[-\log(2)d_{ij}/h]$ with $d_{ij}=0$ handled exactly as one;
a positive ratio too large to represent safely implies an exponential that
underflows to zero. It does not turn the entire requested covariance into its
zero-half-life limit. For long scales, a stable projected representative is

$$
 -\tfrac12 R_q^\top[\,d_{ij}\phi(\log(2)d_{ij}/h)\,]_{ij}R_q,
 \qquad
 \phi(x)=\frac{-\operatorname{expm1}(-x)}{x},\quad\phi(0)=1.
$$

After projection this is algebraically $R_q^\top K(c)R_q/(2c)$; the
constant term vanishes. Evaluate small ratios and $\phi$ stably, avoiding
avoidable overflow, cancellation or a $0/0$ operation. An equivalent stable
kernel is permitted. No finite sentinel replaces either exact endpoint law.
[I, D]

Independent convergence checks summarized in [verification](verification.md)
support retaining ten nodes without a further global increase. Those checks
are examples, not a certificate for all geometries, directions, or confidence
endpoints. The $q\le100$ support ceiling has the same qualification. A
change of quadrature profile must retain positive masses, the uniform-in-half-life
target, common observed/simulated scoring, and an explicit methods-version
record. [D]

## 5. Bounded heuristic confidence-set search

### 5.1 Target and scope

The mathematical finite-simulation inversion is

$$
 \mathcal C_B=\{h\in[0,\infty]:p_B(h)>\alpha\}.
$$

This set may be disconnected. Its hull is a separate summary. Complete
inversion of valid point-null tests has a pointwise confidence-set argument;
the bounded search below is **not a certificate that every component or gap
has been located**. Output from that search is labelled a *detected numerical
confidence set*. Missing an accepted island can remove accepted parameters;
missing a rejection gap can enlarge the reported set. [P, I, D]

The version-0.6 profile retains the following current finite recipe, rather
than asserting an unvalidated simpler or certified alternative. Constants
below are package search policies. They are not specified by the paper and
do not establish a probability-error bound. [D]

### 5.2 Initial grid

Use $t=h/(1+h)$, with $t=1$ representing $h=\infty$.
Construct and sort the union of the following points. [D]

1. The exact endpoints $t=0,1$ and regular points $t=i/256$ for
   $i=1,\ldots,255$.
2. Short-scale points: initialize $h=1/255$; while
   $h>\max(d_{\min}/32,\nu)$, replace $h$ by $2^{-1/8}h$ and
   insert $h/(1+h)$. There are eight probes per doubling in half-life.
   Because the §4 preflight cap of 64 quadrature panels already requires
   $r\ge2^{-63}$ ($d_{\min}\ge2^{-58}$), every geometry passing quadrature
   preflight generates at most **441 short-scale probes**; no separate
   secondary short-probe cap is needed, and excess short-scale work fails
   at §4 preflight before simulation.
3. Long-scale points: define
   $$
   h_{\max}=\min\{10^8,\max\{1024,1/(10^{-4}d_{\min})\}\}.
   $$
   Initialize $h=255\,2^{1/4}$; while $h<h_{\max}$, insert
   $h/(1+h)$ and replace $h$ by $2^{1/4}h$. There are four probes
   per doubling. Compute the bound without avoidable floating-point overflow.

During sorted deduplication, retain a new point $t$ after previous retained
point $s$ only when $t-s>10^{-12}\min(t,1-s)$. Preserve both exact
endpoints, including $t=1$. The finite upper probe bound and the explicit
infinite endpoint are different objects. [D]

### 5.3 Midpoint refinement

Evaluate rank counts on the initial grid. Apply at most **three** refinement
rounds. In each round, operate on the currently sorted grid, using its current
point counts and acceptance states, and collect all proposed midpoints before
evaluating them as a batch. [D]

For adjacent points $t_i<t_{i+1}$, define
$e_i=\min(t_i,1-t_{i+1})$ and the minimum eligible width

$$
 w_i=
 \begin{cases}
 10^{-4}e_i,&e_i>0,\\
 10^{-12}\max\{\min(t_{i+1},1-t_i),\nu\},&e_i=0.
 \end{cases}
$$

Skip midpoint refinement when the interval width is at most $w_i$, or
when its endpoint acceptance states differ. For same-state endpoints, insert
the midpoint if either applicable condition holds: [D]

| Endpoint state | Level condition | Additional local-shape condition |
| --- | --- | --- |
| Both accepted | $\min(p_i,p_{i+1})\le\max(6\alpha,\alpha+0.05)$ | Both outer neighbors exist, $p_i<p_{i-1}$, and $p_{i+1}<p_{i+2}$ |
| Both rejected | $\max(p_i,p_{i+1})\ge0.2\alpha$ | Both outer neighbors exist, $p_i>p_{i-1}$, and $p_{i+1}>p_{i+2}$ |

Merge evaluated midpoints with existing points and apply the same sorted
deduplication rule. Stop early if a round proposes no points. These thresholds
are broad search triggers; in particular, the accepted-cell threshold at 95%
confidence is 0.30. They must not be described as a tight statistical or
numerical definition of “borderline.” [D]

### 5.4 Transition brackets, endpoint bounds, and finite budgets

After midpoint rounds, identify every adjacent pair with different acceptance
states. Bisect each such bracket until its width in $t$ is at most

$$
 \tau_t=\min(10^{-7},10^{-5}\Delta t_{\mathrm{initial}}),
$$

where $\Delta t_{\mathrm{initial}}$ is the width of that transition
bracket at the start of bisection. A midpoint that has already been evaluated
reuses its cached count. Retain the acceptance state at both final bracket
ends. An unresolved numerical point, midpoint stagnation, or unmet accuracy
condition produces incomplete-search behavior below. [D]

Construct components from the located acceptance changes. Round a lower
component endpoint down to the rejected-side lower bracket bound and an upper
endpoint up to the rejected-side upper bracket bound, thereby including the
accepted side of each located transition. If zero is accepted, the first
component starts at zero. Record exact endpoint-test membership separately
from these outer numerical bounds. In particular, `exact_inf_accepted`
records acceptance of the exact $h=\infty$ experiment, whereas
`upper_bound_is_inf` records an infinite numeric enclosure. The two fields
need not agree. An empty detected set is reported explicitly. The hull, when
requested, is the smallest enclosing interval of the detected components; it
is not a replacement for those components or their endpoint-membership flags.
[D]

For every finite transition return both bounds in $t$, both bounds mapped
to normalized $h=t/(1-t)$, and both bounds in requested output units.
Since $dh/dt=(1+h)^2$, $\tau_t$ is not a uniform absolute or relative
precision claim in half-life units. An upper enclosure may be infinite when
the rejected side is **exactly $t=1$**, or when conversion of a finite
$t<1$ overflows in the requested output units. Record which occurred. If
the exact infinite-half-life test rejected, do not print an interval closed
at infinity or describe the unit-root point as accepted: show the infinite
outer enclosure with an excluded infinite endpoint and an endpoint-resolution
note. This limitation concerns the located boundary; it is not evidence that
all sufficiently large finite half-lives are accepted. The same separation
applies if an outer lower bracket reaches zero but the exact zero test
rejected: display an open zero endpoint and preserve that rejected endpoint
test. An infinite numeric bound or a zero outer bound alone never determines
endpoint-test membership. [I, D]

The quadrature panel cap, long-tail cap, deduplication tolerance, three
midpoint rounds, and transition precision define this profile's finite search
budget. There is no additional universal evaluation-count cutoff in the
method definition. Implementations may reject an operation at resource
preflight or impose an explicitly reported resource limit. Resource rejection,
interruption, or failure to meet a required point/bracket accuracy condition
must not be returned as an ordinary finalized confidence set. [D]

Completing all prescribed evaluations yields status `heuristic_complete`,
meaning completion **under this profile**, not completeness of the continuous
set. If the search is interrupted, hits a resource limit, or fails an accuracy
or bracket check, fail with an actionable error by default. A host MAY
optionally return completed point-null evaluations as partial diagnostics with
status `incomplete` and a reason, in which case final confidence components,
hull endpoints, and any claim of nominal confidence-set completion must be
absent (and partial brackets may be exposed only as unresolved diagnostics).
No resumable search state is required; callers may always inspect specific
values with the stateless point-null operation and explicit draws. [D]

## 6. Half-life result and reuse contract

A completed result must include the following, with conventions shared across
host languages. These are semantic requirements; adapters may use idiomatic
container types without changing their meaning. [D]

| Result information | Required meaning |
| --- | --- |
| Sample and geometry | Observation identities, metric, physical units, sample diameter, $q$, and validated positive-rank/block information; a verified lower rank bound with unknown total rank is permitted after a partial eigensolve |
| Simulation and score | Methods/profile version, requested level, actual draw use and identity metadata, and quadrature policy and actual node count where computed |
| Set | Every detected component, optional separate hull, explicit empty/full state, numeric outer bounds and endpoint membership, and whether the set comes from heuristic inversion or a conservative full-domain fallback |
| Endpoint experiments | Evaluations at exact zero and infinity in a completed heuristic search, or an explicit bypass reason for an analytic/policy full-domain result; `exact_inf_accepted` and `upper_bound_is_inf` are separate, with the corresponding zero-endpoint distinction and any resolution limit |
| Full evaluation mesh | **Every distinct point test actually evaluated**, including initial, tail, midpoint and bisection points: $t$, normalized and requested-unit half-life, actual count/draw count if simulated, p-value and its source, acceptance, stage/round, and evaluation status |
| Transition brackets | Initial and final bracket endpoints, endpoint acceptance states, stopping tolerance, mapped half-life bounds, and outward-rounding information |
| Search diagnostics | `heuristic_complete`, structural uninformative, `uninformative_within_profile`, unattainable-rank-level full domain, or `incomplete`; number of rounds/evaluations, applicable limits, and any failure reason |

The exposed mesh must be sufficient to audit the actual reported components.
Returning only the initial grid under an unqualified “grid/p-values” label is
not conforming. All evaluated points remain in the diagnostic trace even if
deduplication excludes a point from the active search mesh. A shortcut result
records only actual work, identifies unevaluated endpoint experiments and their
bypass reasons, and may have an empty trace. Do not fabricate counts or perform
otherwise unnecessary simulations, quadrature construction, or a full
eigendecomposition merely to populate metadata. [D]

Within a prepared operation, reuse normalized geometry, the retained basis,
projected mixture covariances and factors, observed mixture scores, and
point-null evaluations for identical $h_0$, draws and numerical profile.
Evaluate parameter batches and simulation tiles without materializing all
null-by-alternative-by-draw scores. Cache identity must respect each quantity's
dependencies: geometry and covariance factors may be shared across responses
when sample, metric, retained space and profile agree; observed scores and
point counts additionally depend on the response, and simulated counts depend
on draw identity. This requirement does not imply a persistent global cache. [D]

## 7. Spatial transformations

### 7.1 Shared input and output rules

Apply a transformation to each selected numeric response or design column on
one declared analysis sample, using identical geometry and neighbor/cluster
membership across those columns. Preserve fitted observation order and stable
row identity. Define missingness before constructing geometry; do not change
the sample separately for different transformed columns in one regression.
The standalone formula operation follows the common complete-case contract
and restores excluded rows as explicit missing output values. [D]

If a formula is accepted, construct the original model design using its host's
normal contrast and transformation state **before removing an intercept
column**. Standalone and integrated workflows must agree for the same formula,
sample, fitted transformation state and numeric columns. Formula evaluation is
host work; the numerical transformation does not reinterpret formula syntax.
Offsets and unsupported design forms fail explicitly. A standalone
transformation may return rank-reduced numeric columns. The additional
identifiability and separable-constant requirements for fitting the combined
regression are specified in §8; they do not prohibit those standalone outputs.
[D]

All transformations below annihilate a constant input. Compute means and
differences using the stable centering/accumulation rules in
[numerics](numerics.md). Addition of a representable constant and nonzero
rescaling of a numeric input must have the corresponding contrast behavior
within the declared numerical tolerance. Overflow or a non-finite transformed
contrast produces a clear error; it is not silently clipped. [I, D]

### 7.2 Centered LBM-GLS

Using $B=-MdM/2$, let $V$ contain the retained orthonormal eigenvectors
and $\Lambda$ their eigenvalues. The transformation is

$$
 y^*=M V\Lambda^{-1/2}V^\top M y.
$$

This is the symmetric Brownian inverse square root on its retained positive
space, with zero action on discarded directions. The implementation centers
before whitening and recenters the result. Retain eigenvalues **strictly
greater than $10^{-10}$** in the unit-diameter covariance. This absolute
threshold is inherited numerical regularization, not a theorem threshold and
not a relative cutoff. It differs from the positive-rank threshold for
low-frequency diagnostic bases. Report the retained whitening rank and reject
an empty whitening space. Apply the common eigensolver accuracy and ambiguity
rules rather than making an arbitrary basis choice at an unresolved cutoff.
[P, A, D]

The mathematical transformation is invariant to rotations within a wholly
retained eigenspace. Do not construct and store a dense whitening operator
when the same action can be performed from retained factors. A combined
workflow should reuse a suitable full Brownian decomposition for both
whitening and its leading-direction diagnostics. [I, D]

### 7.3 Nearest-positive-distance differences

For observation $i$, let $\delta_i$ be the smallest **strictly positive**
distance $D_{ij}$ to a distinct location. Let $N_i$ contain every observation at
locations tied at that distance under the numerical rule below. Define

$$
 y_i^*=y_i-\frac{1}{|N_i|}\sum_{j\in N_i}y_j.
$$

This extends the paper's nearest-neighbor difference by averaging ties and
explicitly excluding all observations at the query's own location from the
neighbor set. Observations at a selected neighboring location all receive
equal observation weight; do not silently replace them by a single equally
weighted location. Fewer than two distinct locations is an eligibility error.
These co-location and tie choices are package decisions. [P, D]

A distinct location $j$ is included in the numerical tie set when

$$
 D_{ij}-\delta_i
 \le10^{-11}\delta_i.
$$

Evaluate the comparison without avoidable overflow or loss of a positive
distance. This pure distance-relative tolerance is invariant under exact
Euclidean coordinate translations (and requires no coordinate-origin ULP
calculation), while remaining wide enough to preserve ties on ordinary decimal
geographic grids (such as $0.1^\circ$ and $0.01^\circ$ rasters whose binary64
degree differences vary by $\sim 10^{-13}$–$10^{-12}$ relative to the grid
step) and narrow enough to separate distinct near-neighbors. Supplied
coordinates define the geometry; intended ties lost when source coordinates
suffered larger precision loss (for example, rotating unit-spaced coordinates
around an origin of magnitude $10^6$) are not reconstructed silently. [D]

Apply this tolerance only to the nearest-distance comparison. It must not
merge distinct coordinates into co-location groups or change the strict ISO
radius rule. Neighbor indexing/pruning may be conservative, but the final
membership calculation must satisfy the same rule as a direct search. [D]

### 7.4 Isotropic-radius differences

For a positive finite radius $b$, define
$N_i(b)=\{j:j\ne i,\ D_{ij}<b\}$ and

$$
 y_i^*=\begin{cases}
 y_i-|N_i(b)|^{-1}\sum_{j\in N_i(b)}y_j,&|N_i(b)|>0,\\
 0,&|N_i(b)|=0.
 \end{cases}
$$

The radius comparison is **strict**: a pair exactly $b$ apart is excluded.
Other observations at the same coordinates are included for every positive
radius. The isolated-row zero is the specified transform; do not substitute
the untransformed value or an undocumented nearest neighbor. In a subsequent
regression an all-zero transformed row contributes no cross-products.
Geographic radii are in meters and Euclidean radii are in the supplied
coordinate units. Report isolated-row identities/counts. [P, D]

Use an exact metric index or equivalent calculation and accumulate neighbor
contrasts without storing the whole neighborhood graph when possible. A large
radius may legitimately require dense-neighborhood work; the specification
does not promise subquadratic time for every radius. [D]

### 7.5 Cluster demeaning

For a supplied partition $g(i)$, define

$$
 y_i^*=y_i-\frac{1}{n_{g(i)}}\sum_{j:g(j)=g(i)}y_j.
$$

This transformation depends only on group membership and uses equal
observation weights. A singleton becomes zero. A one-observation sample is
permitted for this transformation, although it cannot support the spatial
diagnostics. Missing or invalid cluster labels fail the shared input contract.
No coordinate requirement or within-cluster coordinate constancy is imposed
by this standalone transformation; SCPC's declared-site aggregation has a
different contract. [P, D]

## 8. Combined regression workflow and interpretation

The supported combined workflow requires a levels OLS design with **one
explicit, separable nonzero constant column**, of any name or scale. It need
not be named `Intercept`. Fit that declared levels design and compute the
applicable diagnostics, remove the constant column separately from the slope
columns, transform the fitted response and remaining design with the same
LBM-GLS operator, and fit the transformed regression **without an added
intercept**. The transformation sends the original constant to zero; adding a
new constant after transformation is a different regression, not recovery of
the original intercept. [I, D]

Require the remaining transformed slopes to be jointly identifiable under the
common fitted-design rank rule, preserving their coefficient identities. An
unidentified transformed slope is a clear error, not an arbitrarily dropped
coefficient. A full-indicator or other constant-spanning parameterization
without a separately removable constant column is unsupported by this combined
workflow: explain that the user should refit with an intercept and explicit
contrasts. Do not silently reparameterize the coefficients. Residual SPUR
diagnostics retain their broader constant-span support, and standalone
transformations may return rank-reduced columns. [D]

The levels and transformed fits use the same declared analysis observations
unless a separately named workflow explicitly specifies otherwise. Preserve
coefficient identities, formula contrasts, and fitted transformation state.
Report the transformation, whitening rank, retained diagnostics dimension,
sample identities, and inference profile. Reuse geometry and eligible
decompositions while recomputing model-dependent scores and projections for
each distinct fitted regression. [D]

Müller and Watson (2024), §5, supplies simulation evidence and applied guidance
for transformed regressions; it does not prove validity for every spatial
process or regressor design. The published practitioner decision tree has
recommendation provenance, as recorded in [sources](sources.md), and is not a
general post-selection coverage theorem. A prespecified transformation has its
own inference assumptions. Choosing among analyses using diagnostics on the
same observations adds a selection issue. The combined workflow must not
describe the selected result as having a uniform nominal coverage guarantee
without a separate argument. [R, D]

Return diagnostics and prespecified analysis results with their separate
interpretations. Do not silently turn a diagnostic rejection into a different
estimand, method, or confidence claim. Marginal coefficient intervals are not
simultaneous intervals, and a successful matrix calculation does not validate
the scientific regression assumptions. [D]

## 9. Verification obligations and historical comparisons

Conformance checks must distinguish the point-null score, its finite-rank
calibration, and numerical set search. Verify density-mixture convergence on
independent ordinary and multiscale geometries; endpoint laws including exact
co-location; integer-rank boundary decisions; disconnected/empty/unbounded
output semantics; full evaluation-trace reproducibility; incomplete-search
behavior; and constant, response-scale, distance-unit, row-order and
within-subspace rotation invariances with the distinctions in
[the common reproducibility contract](README.md#7-result-contract-and-reproducibility).
Mathematical operators and retained projectors have their stated deterministic
invariances; angular inference has distributional invariance under a change
of retained basis. Fixed-draw comparison requires the same basis and coloring
convention, or the corresponding rotation of the **colored simulated
realizations**. Reusing identical white-normal arrays with a newly oriented
basis and its new Cholesky factor does not provide that coupling. An
implementation may use a deterministic basis convention without a shared
bespoke anchoring algorithm. Indistinguishable duplicate-row
symmetries are checked through projectors and inference, without inventing
row identifiers or promising bitwise equality. Include actual likelihood
fixtures in addition to synthetic search-function tests. [D]

For transformations, verify direct versus indexed neighbor membership,
nearest-distance ties, exact radius boundaries, co-location, isolated and
singleton rows, whitening rank/cutoff behavior, formula-contrast consistency,
and the transformed regression's absent intercept. A historical-output fixture
and an independent mathematical check serve different purposes. [D]

Keep the author replication's $0.999$ basis proxy, equal-weight finite
mixture, null grid and finite endpoint sentinel as provenance fixtures where
useful. They are not conforming substitutes for this chapter's target and do
not require production compatibility modes. Intentional result changes must
be recorded with their mathematical or numerical reason and applicable
verification evidence. [A, D]
