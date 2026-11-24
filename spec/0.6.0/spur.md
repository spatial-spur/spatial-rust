# SPUR persistence diagnostics

This chapter is normative for methods version **0.6.0**. It defines the adopted
target; the implementation pinned in the [version overview](README.md) reports
0.5.0 until it passes the new conformance requirements. This chapter does not
assert that the baseline code already conforms.

Provenance tags follow [the overview](README.md): **[P]** published definition,
**[R]** paper-recommended tuning, **[A]** author-software convention, **[D]**
package decision, and **[I]** algebraic consequence. MW2024 below denotes
Müller and Watson, *Spatial Unit Roots and Spurious Regression*, identified in
[the source register](sources.md). The shared [numerical profile](numerics.md)
governs geometry, rank and gap tolerances, arithmetic accuracy, resource limits,
and numerical failure reporting. Half-life inference and transformations are
defined separately in [half-life and transformations](half-life-and-transforms.md).

## 1. Supported questions and evidence

**[P, D]** This profile supplies four scalar diagnostics:

| Diagnostic | Response used | Canonical null | Alternative used to construct the statistic |
| --- | --- | --- | --- |
| I(1), or LFUR | A response after removing its constant | Lévy–Brownian covariance shape | One tuned exponential covariance shape |
| I(0), or LFST | A response after removing its constant | Exponential covariance shapes with average correlation at most 0.03 | A tuned Brownian addition to a fixed reference exponential shape |
| Residual I(1) | Response contrasts orthogonal to an eligible OLS design | Conditional Lévy–Brownian disturbance shape | One tuned exponential disturbance shape |
| Residual I(0) | Response contrasts orthogonal to an eligible OLS design | Conditional exponential disturbance shapes with average correlation at most 0.03 | A tuned Brownian addition to the reference disturbance shape |

**[P]** The canonical experiment is Gaussian, with an unspecified positive
variance scale. The definitions below are based on MW2024, equations (16)–(19)
and §4.5. The paper's broader asymptotic interpretations require its sampling,
process and low-frequency conditions; eligibility under a numerical rank check
does not establish those conditions.

**[P, D]** Residual calibration requires the asserted conditional covariance law
given the design. Independence of the disturbance and regressors supplies the
paper's canonical construction; its asymptotic relaxation is described in
Remark 4.5. OLS normal equations, full column rank, or a conditional mean
restriction alone do not establish that covariance law. The implementation
does not claim to test this scientific assumption from fitted residuals.

**[D]** Residual inputs are unweighted, single-response OLS problems with an
eligible numeric design spanning the constant. They may be evaluated in batches
as specified in §8. Weighted, offset, penalized, nonlinear and IV residual
diagnostics are outside this profile. Their mere possession of coefficients and
residuals does not make them supported inputs. An adapter must apply the
estimator and sample contracts in [the overview](README.md) before calibration.

**[P, D]** Rejecting one persistence null does not prove the other null. Both
tests can reject, and neither can reject. The outputs are test results under
their respective assumptions, not an exhaustive classification into two models.
Selecting a subsequent regression or transformation using these outcomes does
not receive a uniform post-selection coverage guarantee from this chapter.

## 2. Common geometry, projections, and dimensions

### 2.1 Distances and contrast spaces

**[P, D]** Let there be $n$ analysis observations on the common retained sample.
Use the declared metric and exact co-location rule in [the overview](README.md).
Let $D$ be the pairwise distance matrix divided by the sample diameter, so
$0\le D_{ij}\le1$. The diameter must be positive and finite. For SPUR
diagnostics, repeated observations are retained unless the caller has explicitly
prepared another analysis unit; no implicit averaging, jitter or nugget is
introduced.

Define

$$
K(c)_{ij}=\exp(-cD_{ij}),\qquad
K(\infty)_{ij}=\mathbf 1\{D_{ij}=0\},
$$

and

$$
f=\frac{1}{n(n-1)}\sum_{i\ne j}\mathbf 1\{D_{ij}=0\},\qquad
\bar\rho(c)=\frac{1}{n(n-1)}\sum_{i\ne j}\exp(-cD_{ij}).
$$

**[I]** Distinct observations at identical coordinates have unit covariance in
this kernel at every decay, including the infinite-decay limit. Declaring this
benchmark does not assert that observed replicate measurements cannot contain
idiosyncratic noise; that would be a separate model component.

**[P, I]** For levels tests, use
$M=M_1=I-\mathbf1\mathbf1^\top/n$. For residual tests with a full-rank design
$X$, use $M=M_X=I-X(X^\top X)^{-1}X^\top$. The formula defines the
orthogonal projector; its implementation must use the stable factorization
required by [the numerical profile](numerics.md), rather than forming unstable
normal-equation inverses.

**[D, I]** Residual designs must span the constant within the shared numerical
criterion. This support restriction removes the origin-dependent terms of the
Brownian covariance and makes the test invariant to adding a constant. It does
not declare every possible no-intercept residual-test formulation invalid.
The design must include the evaluated fitted columns and transformation state,
not columns reconstructed under different knots, contrasts or sample selection.

### 2.2 Brownian basis and projected covariances

**[P, I]** Define the projected Brownian matrix

$$
B_X=-\tfrac12 MDM.
$$

Let $R\in\mathbb R^{n\times q}$ have unit-norm orthogonal columns spanning its
leading $q$ eligible positive directions. In exact arithmetic
$MR=R$, $R^\top R=I_q$, and $R^\top\mathbf1=0$. Define

$$
z=R^\top My=R^\top y,\qquad
\Omega_B=R^\top B_XR,\qquad
\Omega_K(c)=R^\top MK(c)MR.
$$

For every eligible projected covariance define its unit-trace shape
$H(\Omega)=\Omega/\operatorname{tr}(\Omega)$, and write
$H_B=H(\Omega_B)$ and $H_K(c)=H(\Omega_K(c))$. Compute positive traces and
normalize safely. These shapes remove the arbitrary covariance scale from
the numerical representation as well as from the statistical model. **[D, I]**

The paper scales its eigenvector columns by $\sqrt n$. Consistently using unit
columns rescales the observation and covariance coordinates without changing
the ratio or inference. Implementations must compute the covariance of the
actual retained directions; an approximate numerical tie is not permission to
replace that covariance by an unrelated diagonal matrix.

**[D]** Complete numerical eigenvalue blocks are retained. A requested cutoff
that splits a block is rejected with the requested $q$ and an actionable
description of the eligible dimension. Rank thresholds, gap thresholds and
the basis policy are fixed by [the numerical profile](numerics.md).
A deterministic implementation basis serves reproducibility within its
recorded numerical profile. An orthogonal change of basis must transform the
observation and every covariance consistently. The profile does not require
projector-anchor machinery solely to reproduce the same realized finite-draw
p-value under every row permutation or linear-algebra provider.

### 2.3 Dimension policy and nonzero responses

**[R, D]** The default is $q=15$. Require integer $2\le q<n$, sufficient
positive Brownian rank, and, for residual tests,
$q\le n-\operatorname{rank}(X)$. The requested dimension is not silently
reduced to fit a small sample or a split eigenvalue block. Failure reports the
constraint and the available dimension information.

**[I, D]** $q=1$ has no scale-free covariance-shape information. A larger
eligible $q$ defines a canonical calculation but does not automatically
inherit a fixed-low-frequency asymptotic approximation when it approaches full
rank. $q$ must be prespecified or selected by an outcome-independent
geometry/design rule. Outcome-driven searches over $q$ are not covered by
the reported marginal calibration.

**[D]** Emit one higher-frequency advisory per operation if $q\ge40$ or all
numerically positive directions are used. Record both conditions separately.
These are interpretation advisories, not theorem boundaries.

**[D]** A nonfinite response, a zero response contrast $z$, or a contrast
indistinguishable from zero under the shared relative numerical test is an
error. Report the affected response and reason. Do not assign a statistic by
dividing two zeros, and do not confuse this error with the uninformative
covariance-shape result in §7. Centering, scaling and residual projection must
preserve the invariances and arithmetic safeguards in [the numerical profile](numerics.md).

## 3. I(1) statistic and parameterization

**[D, I]** Parameterize alternative strength by $w\in[0,1]$. In the
unit-diameter geometry, define the contrast-scale reference decay directly from
the already-constructed projected Brownian matrix $\Omega_B=R^\top B_X R$:

$$
c_{\mathrm{ref}}=\frac{n}{2\operatorname{tr}(\Omega_B)},
\qquad
c(w)=c_{\mathrm{ref}}\,\frac{w}{1-w}\quad(0<w<1),
\qquad
H_{1,a}(w)=H_K(c(w)).
$$

This is a direct compact parameterization of every positive finite decay,
requiring no average-correlation root solve at any tuning candidate. Because
$\operatorname{tr}(\Omega_B)/n$ is invariant under orthogonal basis changes and
uniform site duplication, and scales linearly with contrast-subspace distance
(including residual geometries where $X$ absorbs distant outliers), $c_{\mathrm{ref}}$
keeps the 50%-power decay well inside the interior of $(0,1)$ across arbitrary
spatial scales. Define the exact endpoints $H_{1,a}(0)=H_B$ and
$H_{1,a}(1)=H_K(\infty)$. The weak endpoint is the projected exponential
shape limit as $c\downarrow0$; it is not the zero contrast covariance
$M\mathbf1\mathbf1^\top M$. The strong endpoint uses the exact co-location
covariance, with no finite decay sentinel or hardware-underflow stopping rule.
For a selected $0<w_a\le1$, define the reported statistic

$$
T_1(z;w_a)
=\frac{z^\top H_B^{-1}z}
       {z^\top H_{1,a}(w_a)^{-1}z}.
$$

**[P, D]** Reject for large values. Its canonical null is
$z\sim N(0,\sigma^2 H_B)$, $\sigma^2>0$. The alternative used for tuning
is $N(0,\sigma^2 H_{1,a}(w_a))$. Because variance scale is unrestricted,
normalizing these covariances does not change the null or alternative
experiments. Residual I(1) uses the same definitions with the residual $M$,
$R$ and $z$.

**[A, D, I]** Unit-trace ratios replace the author software's scale-dependent
reported ratio. For a fixed finite alternative decay $c_a=c(w_a)$,

$$
T_1
=\frac{\operatorname{tr}(\Omega_B)}
        {\operatorname{tr}(\Omega_K(c_a))}\,T_{1,\mathrm{author}}
=\frac{2c_a\operatorname{tr}(\Omega_B)}
        {\operatorname{tr}(\Omega_K(c_a))}\,T_{1,\mathrm{paper}}.
$$

The paper expression here uses
$\Omega_{\mathrm{LTU}}(c)=\Omega_K(c)/(2c)$ and the same Brownian
normalization and distances. These are positive, response-independent factors;
critical values scale identically, preserving ordering, p-values and decisions
for the same selected alternative and draws. The exact endpoint definitions
stand on their own; no finite-rate conversion is asserted at $c=\infty$.
There is no separate legacy-normalization production mode.

**[D]** Report the selected `alternative_strength` $w_a$ with meaning
`scaled_compact_decay`, its contrast-scale anchor $c_{\mathrm{ref}}$, its
achieved average correlation, and $c_a$ with its units, using an exact-limit
marker when $w_a=1$. A finite $c_a$ multiplies diameter-normalized distance;
for absolute diameter $L$, its decay per absolute distance unit is $c_a/L$.
These are statistic-tuning parameters, not estimates of the observed response's
persistence. Selection depends on geometry, the eligible design and tuning
draws, not on $z$.

## 4. I(0) statistic and composite null

### 4.1 Reference and alternative

**[R, D]** The canonical null bound is $\rho_{\max}=0.03$, and the reference
correlation for choosing the statistic is $\rho_\star=0.001$. If
$f<\rho_\star$, solve $\bar\rho(c_\star)=\rho_\star$. Otherwise set
$c_\star=\infty$. The latter is the adopted attainable-limit extension;
the solver must not claim to have found a finite rate for an unattainable
correlation.

Write $H_\star=H_K(c_\star)$. For strength $w\in[0,1]$, define
$H_{0,a}(w)=(1-w)H_\star+wH_B$. For a selected $0<w_a\le1$, report

$$
T_0(z;w_a)
=\frac{z^\top H_\star^{-1}z}
       {z^\top H_{0,a}(w_a)^{-1}z}.
$$

**[P, D]** The Brownian mixture supplies the alternative shape. The weak
endpoint $w=0$ equals the reference, and $w=1$ is exactly Brownian. The null
used to tune $w_a$ is the reference $H_\star$; the null family used for final
inference is wider. The same $T_0$, including its fixed reference and selected
alternative, must be evaluated under every null covariance. Do not retune the
statistic at each nuisance point or for each observed response.

**[A, D, I]** This bounded strength replaces the author software's unbounded
covariance multiplier. Writing $\Omega_\star=\Omega_K(c_\star)$, the
matching author alternative for $0<w_a<1$ is
$\Omega_\star+\gamma_a\Omega_B$, where

$$
\gamma_a=\frac{w_a}{1-w_a}
          \frac{\operatorname{tr}(\Omega_\star)}
               {\operatorname{tr}(\Omega_B)},\qquad
T_0=(1-w_a)T_{0,\mathrm{author}}.
$$

For finite $c_\star$, the paper's Brownian amplitude satisfies
$\gamma_a=2c_\star g_a^2$ under its LTU normalization. These maps explain
provenance; implementations use the bounded mixture directly, without forming
an overflowing $\gamma_a$ near $w=1$. The positive factor preserves statistic
ordering for the same alternative. The exact Brownian endpoint needs no finite
multiplier. Report `alternative_strength` $w_a$ with meaning
`brownian_mixture_share`; no duplicate author-parameter field or compatibility
calculation is required.

### 4.2 Continuous target and adopted numerical family

**[P, D]** If $f<0.03$, solve $\bar\rho(c_{0.03})=0.03$. The continuous
canonical null is

$$
\mathcal C_0=\{c:c\ge c_{0.03}\}\cup\{\infty\},\qquad
z\sim N(0,\sigma^2\Omega_K(c)).
$$

If $f=0.03$, the null consists only of $c=\infty$. If $f>0.03$, the
specified null is infeasible and the operation errors with the attained floor
and number of distinct locations. This is not a rejection of the I(0) null by
a statistical test.

**[P]** The population p-value target is

$$
p_0(z)=\sup_{c\in\mathcal C_0}
 P_c\{T_0(Z)\ge T_0(z)\}.
$$

**[D]** Version 0.6.0's default computation uses the shared bounded decay mesh
in [the numerical profile](numerics.md) at $\rho_0=0.03$: the 1.05 interior
decay mesh from $c_{0.03}$ to $c_{\max}$, its 1.20 finite decay-tail probes up
to $\max(c_{\max},36/d_{\min})$, and the exact infinite-decay limit $K(\infty)$.
Exceeding the declared interior or tail budget produces an explicit preflight
error rather than silently coarsening the mesh.

For the evaluated mesh $\mathcal G_0$, return

$$
\widehat p_{0,\mathcal G}(z)=\max_{c\in\mathcal G_0}\widehat p_c(z),
\qquad
\widehat k_{0,\alpha}=\max_{c\in\mathcal G_0}\widehat k_{c,\alpha}.
$$

**[D]** Record `calibration_status = finite_grid` and
`probability_method = monte_carlo`. It is a finite-mesh
approximation to the continuous target, not a continuous-size certificate.
The point-null Monte Carlo argument in §6 covers a null covariance included in
the mesh; it does not bound omissions between mesh points or beyond a finite
tail endpoint. Neither small entrywise correlations nor a stated grid ratio
alone establishes a projected-tail-probability bound.

## 5. Reference-experiment tuning

### 5.1 What is tuned

**[R, D]** Choose an alternative targeting 50% power for a 5% reference-null
test. Use one independent tuning block of $B_{\rm tune}$ standard Gaussian
directions, with default $B_{\rm tune}=100{,}000$. Its count may differ from
the inference count $B$ in §6; equal counts are not a theoretical requirement.
Require $B_{\rm tune}\ge19$, so the 5% reference rejection threshold is
attainable. A smaller requested tuning block is an input error explaining that
minimum. Do not use the observed response to choose the alternative.

For reference shape $H_0$, candidate $H_a$, and common tuning directions
$e_b$, form null draws $L_0e_b$ and alternative draws $L_ae_b$, where
$L_0L_0^\top=H_0$ and $L_aL_a^\top=H_a$. For these simulations, use the
positive-diagonal lower-triangular Cholesky factors in the implementation's
retained basis.
Equivalent algebra may reuse factors and quadratic forms if it preserves the
specified colored vectors for supplied draws. Evaluate the candidate's same
quadratic-form ratio on both sets. The tuning cutoff uses the same rank rule
as final inference in §6: the $k$-th largest null statistic, with
$k=\lfloor(B_{\rm tune}+1)/20\rfloor$, calculated by checked integer
division. The estimated reference power is the fraction of alternative
statistics **strictly above** that cutoff.

**[D]** This reference calculation chooses the statistic; it does not supply
the final inference draws or replace the final composite-null calibration.
One rank-cutoff implementation is sufficient for both operations. The author
software's interpolated type-7 tuning quantile is not retained as a separate
production convention.

**[D]** For I(1), $H_0=H_B$. For I(0),
$H_0=H_\star$, even though final I(0) inference protects the
larger nuisance family. Consequently, the reported I(0) tuning power is not the
power of its final composite-null test.

### 5.2 Bounded strength search

**[D]** Both diagnostics use one search on the bounded strength interval.
The weak limit $w=0$ has an analytically constant statistic and zero strict
exceedance power. It supplies a lower reference but is never selected as an
alternative. The search uses at most **20 interior candidates**, including
$w=1/2$, plus the exact strong endpoint $w=1$:

1. Evaluate $w=1$ and $w=1/2$. Keep every successfully evaluated positive
   strength with its estimated reference power.
2. If any retained power is within 0.01 of 0.5, select the retained candidate
   closest to 0.5; break ties by the smaller strength.
3. Otherwise, in the ordered points $0,1/2,1$, find adjacent endpoints with
   lower power below 0.5 and upper power above 0.5. If no such bracket exists,
   select the retained positive candidate closest to 0.5 and report
   `target_unattained`, with reason `no_evaluated_bracket`. This is not a proof
   that another strength or the population power could not attain the target.
4. With a bracket, evaluate its arithmetic midpoint. Stop if a retained
   candidate meets the power tolerance. Otherwise replace the lower endpoint
   when the midpoint power is below 0.5, or the upper endpoint when it is above
   0.5. Repeat until the tolerance is met or the 20-interior-candidate budget
   is exhausted. Select the retained positive candidate closest to the target,
   again breaking ties by smaller strength. A missed tolerance has status
   `target_tolerance_not_met`, not a convergence claim.

Record the actual evaluated strengths, powers, selected candidate, candidate
count and stopping reason. The analytic weak endpoint must be identified as
such; it is not an invented simulation evaluation. The selected alternative
must be one of the successfully evaluated positive strengths. The exact strong
endpoint is a permissible selected alternative for either diagnostic.

**[D]** Covariance shapes close to the weak endpoint must use the stable
projected-exponential construction in [the numerical profile](numerics.md),
including `expm1` where needed, before trace normalization. Arithmetic failures,
invalid factors or unresolved geometry remain errors. Failure to reach a
preferred reference power after valid evaluations does not invalidate the
independently calibrated test and is reported as a tuning status.

This is a bounded tuning heuristic, not a proof of monotone power or a global
optimizer. Finite-draw power can be stepwise, and the chosen exponential shape
path need not be assumed to have globally monotone power. The old correlation
anchors, special upper starts, 200-step doubling/halving phases and
underflow-based endpoint detection are not part of this profile.

**[D]** If the reference and candidate shapes agree within the shared shape
tolerance at a trial, assign zero strict-exceedance tuning power for that
trial as a conservative profile convention. Call the statistic algebraically
constant only under exact proportionality. A weak trial or a proportional
strong endpoint alone does not establish that all intermediate I(1) shapes
are uninformative. Apply §7 only to a structurally uninformative I(0) family or
the finally selected pair.

### 5.3 Reporting and warnings

**[D]** The canonical diagnostic name is `reference_tuning_power`. It is the
estimated reference-experiment power for the **selected** alternative. The
0.5.0 field `best_tuning_power` denotes this quantity; a host may document a
name mapping, but 0.6.0 requires neither duplicate fields nor a compatibility
algorithm. The word “best” must not be interpreted as a continuous optimum or
final composite-test power.

Report the selected strength and its diagnostic-specific meaning, target 0.5,
reference size 0.05, achieved estimate, tuning draw count, number of evaluated
candidates and
tuning status. Emit a single low-power advisory per operation when an achieved
estimate is below 0.45, naming the affected diagnostic. Always expose a missed
tolerance or budget status even when the estimate is above that warning
threshold. An uninformative result has its own status and no purported
50%-power parameter.

## 6. Simulation, p-values, and critical values

### 6.1 Draw contract and reproducibility

**[D]** The default is Monte Carlo inference with $B=100{,}000$ columns and an
independent tuning block with $B_{\rm tune}=100{,}000$ columns. Users may
request any positive integer $B$ and any integer $B_{\rm tune}\ge19$ within
the shared resource limits. Defaults and explicit controls for these counts
must be documented independently; changing the inference count need not
change the tuning count. Supplied blocks have $q$ rows but may have different
column counts. The API reports both actual counts; a diagnostic uses
$B+B_{\rm tune}$ generated columns in total.

Each block must have independent isotropic standard Gaussian directions, or
equivalently independent directions uniform on the unit sphere. Both blocks
must be independent of the observed outcomes, and the inference block must be
independent of the tuning block. The host RNG/seed policy is fixed in
[the overview](README.md). Equal seeds in different RNG implementations are not
a cross-host identity guarantee. Identical supplied blocks support numerical
parity comparisons when the retained bases and coloring coordinates are also
aligned. Otherwise orthogonal basis changes preserve the simulation law, not
necessarily the realized Monte Carlo p-value; compare invariant quantities
and assess simulation variation as specified in [verification](verification.md).

**[I, D]** Multiplying each column by any nonzero scalar preserves these
homogeneous ratio calculations. Normalize finite nonzero columns by their
maximum absolute entry for arithmetic stability. All-zero or nonfinite columns
are errors. The arithmetic check cannot certify a user-supplied block's
distribution or independence; supplying arbitrary deterministic columns does
not inherit the probabilistic finite-simulation guarantee.

Generated tuning and inference blocks must use disjoint draws from the
declared generator; do not initialize both blocks at the start of the same
seeded stream. An explicitly documented independent-stream scheme is also
permitted. **[D]**

**[D]** Reuse the inference block across nuisance values and the tuning block
across candidate alternatives. The same blocks may be reused across diagnostics
or response columns sharing a sample. This common-random-number convention does
not require independence between the resulting p-values. It also does not
license reusing inference columns as tuning columns.

### 6.2 Ranks, cutoffs, and equality

**[D, I]** At one null covariance $\Omega$, let $T_b$ be the statistic on
the recolored inference draw and $T_{\mathrm{obs}}$ the statistic on the
observed response. Recolor with the positive-diagonal lower Cholesky factor,
using the same factor convention as tuning. Return

$$
\widehat p_\Omega
=\frac{1+\#\{b:T_b\ge T_{\mathrm{obs}}\}}{B+1}.
$$

Equality is included in the tail. Reject at nominal size $\alpha$ when
$\widehat p_\Omega\le\alpha$. Report cutoffs at 1%, 5% and 10%. With
$k_\alpha=\lfloor\alpha(B+1)\rfloor$, the cutoff is the
$k_\alpha$-th largest simulated statistic if $k_\alpha\ge1$, and
$+\infty$ otherwise. Reject **strictly above** that cutoff. For the default
1%, 5% and 10% levels, calculate the ranks by checked integer division:
$\lfloor(B+1)/100\rfloor$, $\lfloor(B+1)/20\rfloor$, and
$\lfloor(B+1)/10\rfloor$, respectively. Floating-point multiplication by
the displayed decimal must not shift a rank at equality. The package's default
decision flags use these exact rational levels and integer counts; displayed
p-values are their numeric representation, not a separate rounding rule.

**[I]** For the same simulated statistics these definitions agree, including
ties:

$$
\widehat p_\Omega\le\alpha
\quad\Longleftrightarrow\quad
T_{\mathrm{obs}}>\widehat k_{\Omega,\alpha}.
$$

For I(0), maximize both the point-null p-values and the point-null cutoffs over
the **same** nuisance mesh, with the same selected statistic and inference
block. This preserves the corresponding finite-mesh decision equivalence.
Do not compute a p-value and cutoff under different meshes, different
alternative choices or separately redrawn simulations.

**[D]** The minimum attainable p-value is $1/(B+1)$. Report any of the default
nominal levels at which rejection is impossible. Infinity is a genuine
unattainable-threshold value; a host lacking an infinity representation must
return an explicit flag rather than a finite sentinel.

### 6.3 Meaning of the guarantee

**[I, D]** Conditional on eligible geometry/design and independent tuning, the
inclusive rank construction is conservative under the specified point-null
Gaussian experiment when observation and inference draws are exchangeable in
the relevant statistic. This statement includes simulation randomness. It is
not a guarantee for arbitrary supplied deterministic blocks or outcome-driven
selection of the method.

For an I(0) covariance included in the evaluated mesh, taking the maximum
p-value cannot increase rejection relative to its correctly calibrated
point-null p-value. No such argument fills an unexamined part of the continuous
null family. Report Monte Carlo resolution and the numerical calibration status
separately from claims about the scientific model or nuisance-search coverage.

**[I, D]** Rank resolution $1/(B+1)$ is not a bound on Monte Carlo probability
error. Conditional on a fixed observed statistic and selected alternative,
the independent point-null exceedance count has a binomial distribution. A
tail probability near 0.05 with $B=100{,}000$ has standard error about 0.00069,
despite a rank step near 0.00001. Report enough information to recover the
point-null counts and identify the simulation method; do not describe decimal
resolution as probability accuracy. A standard error at one maximizing mesh
point is not an uncertainty bound for the composite maximum or the tuning
search.

## 7. Uninformative covariance shapes and failures

**[I, D]** A nonzero response can provide no discrimination for the selected
scale-invariant experiment if its null and alternative covariance shapes are
proportional, or indistinguishable under the shared shape tolerance. Such a
case is a valid conservative result: return `pvalue = 1`, all reported critical
thresholds $+\infty$, and status `uninformative_within_profile`, with a
reason identifying the shape comparison. Do not manufacture
a tuned 50%-power alternative. If no alternative was selected, report its
parameter, the raw ratio and `reference_tuning_power` as unavailable, with
the reason. If a valid alternative was already selected, report its defined
ratio as uninformative; call the ratio constant only under exact
proportionality, not merely agreement within a numerical tolerance.

For I(0), test the reference shape against the Brownian shape before tuning:
exact proportionality makes every mixture in this specified alternative
family uninformative; agreement within tolerance invokes the declared
conservative profile decision. For I(1), run the bounded search. Matching a few
anchors or the strong endpoint alone is not that proof. Apply the uninformative
rule to the finally selected pair when it is indistinguishable at the declared
tolerance. A conservative result for that selected pair is not a conclusion
about every untested shape.

Validate the response projection, null feasibility and all required arithmetic
conditions before returning an uninformative result. A reusable geometry-only
calibration status must not mask a bad response column or unsupported input.

**[D]** In contrast, zero response projection, insufficient rank, a split
eigenspace, infeasible I(0) correlation range, materially invalid covariance,
an unrepresentable calculation, or exhausted required mesh resources are
errors with the corresponding diagnostic. Do not turn these failures into
`pvalue = 1`, silently reduce $q$, add a ridge/nugget, or omit a requested
null covariance. Details follow [the numerical profile](numerics.md).

## 8. Multiple responses and no-cointegration interpretation

### 8.1 Multiple response columns

**[D]** A matrix $Y=[y_1,\ldots,y_m]$ means a batch of $m$ **marginal
scalar tests**. On a common retained sample and common design, reuse geometry,
basis, projected covariances, tuning and simulation blocks, and calculate each
$z_j=R^\top My_j$ and observed statistic separately. Preserve response
identity and order. Any nonfinite or zero projected response is identified as
an error before returning a purported complete successful batch; no column is
silently omitted.

The common-sample missingness and row-identity rules are those in
[the overview](README.md). Columns requiring different samples or designs are
separate diagnostic problems; shared numerical arrays do not authorize a
different analysis sample for each column without recomputing its geometry.
A host interface restricted to one response may implement the same batch by
explicit scalar calls; it must not claim a joint test merely because an
internal kernel accepts a matrix.

**[D]** Batch output provides no family-wise, false-discovery, joint-vector or
simultaneous-confidence guarantee. It must not report the minimum p-value as a
valid joint p-value. Any multiplicity adjustment is a distinct declared
procedure, outside the default SPUR profile. Dependence between responses does
not by itself invalidate a correctly specified marginal experiment.

### 8.2 Residual I(1) as a no-cointegration diagnostic

**[P]** MW2024 Remark 4.2 supplies a no-cointegration interpretation for residual
I(1) under its joint canonical model. The stochastic columns of $(Y,X)$,
excluding deterministic constants, are full-rank linear combinations of
independent Lévy–Brownian fields. A population linear projection then leaves
a Lévy–Brownian disturbance independent of the regressor fields. Conditioning
on $X$ and projecting away its span produces the residual I(1) experiment
above without needing the unknown population coefficient.

**[P, D]** This interpretation concerns that specified no-cointegration null
and its stated asymptotic extension. It does not make every OLS residual test
a valid no-cointegration test for arbitrary endogenous, heteroskedastic or
misspecified joint processes. A rejection establishes departure from the tested
benchmark; it is not uniquely diagnostic of cointegration under every other
possible model. Non-rejection does not prove absence of cointegration.

**[D]** The profile does not estimate a multivariate cointegrating rank,
automatically search over response normalizations, or construct a joint
confidence region for cointegrating vectors. Applying scalar tests to several
normalizations remains several tests. A confidence set obtained by testing
prespecified candidate linear combinations requires its own declared inversion
procedure and is not implicitly supplied by a batch call.

## 9. Required result meaning and numerical evolution

### 9.1 Result record

**[D]** Every successful diagnostic must identify the methods version, test
kind and levels/residual variant; observed sample and response identities;
requested/used $q$; validated positive-rank and retained-tie information; and
whether all positive directions were used. A partial eigensolve may report a
validated rank lower bound, unknown total rank and unknown all-directions-use
status; metadata must not force a full decomposition solely to fill those
fields. Report the unit-trace statistic normalization,
p-value, critical thresholds, alternative strength and its meaning, the I(1)
decay/correlation when applicable,
`reference_tuning_power`, tuning status, numerical calibration status, draw
counts and reproducibility metadata.

For I(0), expose the actual evaluated nuisance mesh with its decay parameters,
achieved average correlations and exact-limit markers, together with the
point-null p-values and cutoffs. Average correlations alone are insufficient
to identify tail evaluations that round to the co-location floor. The returned
mesh size varies with geometry and co-location floor; do not claim a fixed
point count. Shared execution and accuracy fields follow
[the numerical profile](numerics.md).

### 9.2 Deterministic quadratic-form evaluation

**[I, D]** Monte Carlo is the adopted default in 0.6.0, not a mathematical
necessity of LFUR or LFST. For a fixed statistic
$T(z)=(z^\top A_0z)/(z^\top A_az)$ and $z=L\xi$, its upper tail is the
Gaussian quadratic-form event

$$
\xi^\top L^\top(A_0-tA_a)L\xi\ge0.
$$

A future general Davies/Imhof-type evaluator may implement the same population
target only after the accuracy, failure-reporting, independent-reference,
statistical-equivalence and representative-performance gates in
[verification](verification.md) are satisfied. It must be identified in
execution metadata and versioned if the adopted finite conventions change.
The specialized SCPC probability routine is not presumed to be a general
indefinite-quadratic-form solver. No deterministic rewrite is required for
0.6.0 conformance, and it does not automatically solve continuous nuisance
optimization or the separate half-life mixture statistic.

**[D]** Implementations should reuse immutable geometry, basis and projected
covariances where the problem is identical. Reuse must preserve the distinct
levels and residual contrast spaces and every reported response's sample.
Changing storage, tiling or factorization may preserve this specification;
changing the benchmark, retained subspace, calibration mesh or inferential
scope requires an explicit methods decision.
