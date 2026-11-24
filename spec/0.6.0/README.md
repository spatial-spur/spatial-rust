# Spatial SPUR methods specification 0.6.0

**Status: adopted specification and implementation target.** This version
replaces the mathematical and supported-use descriptions of earlier editions
when an implementation declares conformance to 0.6.0. The reviewed baseline is
commit [c41aaaca50559d2750bf10b2d5e8d2d1aac2746d](https://github.com/spatial-spur/spatial-rust/commit/c41aaaca50559d2750bf10b2d5e8d2d1aac2746d),
which reports methods 0.5.0 and does **not** yet satisfy all requirements below.
Adding this document does not change the runtime method version. The
[implementation gates](verification.md) govern that change.

## 1. Authority and document map

This directory is one specification. The following chapters are incorporated
normatively; no historical specification is needed to define the method:

| Chapter | Authority |
| --- | --- |
| [This document](README.md) | Scope, assumptions, fitted samples, geometry, errors, outputs, and versioning |
| [SCPC](scpc.md) | Coefficient calculations, calibration families, dimension selection, and Nyström |
| [SPUR diagnostics](spur.md) | I(0), I(1), residual tests, tuning, and simulation inference |
| [Half-life and transformations](half-life-and-transforms.md) | Density mixture, inversion, spatial transformations, and combined workflow |
| [Numerical profile](numerics.md) | Finite meshes, accuracy criteria, convergence, reproducibility, and resource policy |
| [Decisions and provenance](decisions.md) | Adopted choices, their origin and rationale, and coverage of the audit |
| [Verification and migration](verification.md) | Required evidence and outstanding baseline implementation changes |
| [Sources](sources.md) | Paper and author-software identifiers |

MUST and MUST NOT are conformance requirements. SHOULD denotes a recommendation
whose departure needs a documented reason. MAY denotes an allowed choice.
Requirements in method chapters define the mathematical object; numerical
details explicitly delegated to numerics.md are defined there. An unresolved
contradiction between chapters blocks conformance; it is not permission to pick
the easier requirement. The decision register is explanatory unless it
explicitly adopts a requirement reflected in a normative chapter.

Provenance tags distinguish:

- **[P]** a published mathematical definition or result, under its assumptions;
- **[R]** a published numerical/design recommendation;
- **[A]** a convention in author software or empirical replication;
- **[D]** a package decision, including its supported numerical profile;
- **[I]** an algebraic consequence of the declared definitions.

A paper recommendation is not a necessary theorem condition. Author software
does not override a published definition. A more accurate implementation of a
paper target need not reproduce a historical approximation. Source citations
identify provenance; they do not import unspecified algorithms into this version.

## 2. Notation

| Symbol | Meaning |
| --- | --- |
| n | Number of fitted observations or observations in a SPUR analysis |
| m | Number of declared SCPC spatial sites |
| A | n by m incidence matrix, one unit entry per observation |
| X, Z | Actual fitted structural regressors and instruments, in fitted row order |
| F | Orthogonal projection off the fitted absorbed-effect space; identity without FE |
| a_j | Bread-adjusted observation direction for coefficient j |
| e_hat | Actual fitted structural residuals |
| q | Number of retained spatial directions |
| D, d_ij | Distance matrix and its entries, in the units stated for the operation |
| Delta | Positive sample diameter in external distance units |
| K(c) | Matrix with entries exp(-c d_ij) |
| f | Fraction of off-diagonal pairs at exactly the same location |
| B | Number of inference draws, and separately the tuning-block length when used |
| alpha | Rejection probability; a confidence level is 1-alpha |

Matrix formulas describe mathematical operations, not a requirement to
materialize dense projectors or invert matrices explicitly.

## 3. Supported calculations and scientific assumptions

### 3.1 Host support

The shared core owns the numerical methods. Host adapters own verified model
extraction, row identity, formula evaluation, and presentation. **[D]**

| Fit or workflow | Python | R | Stata |
| --- | --- | --- | --- |
| Unweighted OLS | Statsmodels OLS and PyFixest feols | lm and fixest feols | regress with verified estimation sample and fit |
| Ordinary unweighted 2SLS | PyFixest Feiv with stored within X and Z | fixest feols with retained structural and first-stage data | ivregress 2sls without absorb |
| Categorical absorbed effects | Supported PyFixest one-way and multiway fits | Supported fixest one-way and multiway fits | Unconditional areg slopes with one plain absorb variable |
| Site clustering | One explicit site label per observation | Same | cluster() or one fitted cluster variable with constant coordinates within site |
| SPUR residual tests | Unweighted OLS design spanning a constant | Same | Numeric OLS regressors with an intercept |
| Combined transformed-regression workflow | OLS with a separable constant column and identifiable transformed slopes | Same | Numeric OLS regressors with an intercept and identifiable transformed slopes |
| Standalone transformations | Evaluated numeric variables/formulas supported by the host | Same | Supported numeric variables |
| Weighted, penalized, robust-loss, GLM, nonlinear fits; offsets | Unsupported | Unsupported | Unsupported |

Adapters MUST identify the fitted estimator, not merely accept an object whose
class resembles OLS. Unsupported estimators MUST be rejected before inference,
including regularized fits and nontrivial weights. An adapter MUST NOT refit an
unsupported estimator as unweighted OLS and present the result as inference for
the supplied fit. Formula offsets are unsupported even if they can be evaluated.
Unsupported FE encodings, varying slopes, incomplete first-stage state, or
ambiguous estimation samples require descriptive errors.

This matrix specifies executable support. Scientific assumptions remain the
user's responsibility; numerical validation cannot establish them.

### 3.2 Evidence and guarantee levels

| Workflow | Interpretation that MAY be claimed | Conditions and limits that MUST accompany it |
| --- | --- | --- |
| SCPC Gaussian benchmark | Calibration for the specified Gaussian ratio family at evaluated covariances | Actual family, numerical status, and finite-mesh limitation must be reported |
| Unconditional regression SCPC | Relevant paper asymptotic justification | Identified coefficient, valid score representation, and applicable spatial weak-dependence/design conditions |
| Conditional OLS SCPC | Conditional Gaussian benchmark protection together with location-family protection | The stated sign/site covariance model; not arbitrary conditional heteroskedastic or within-site covariance |
| Conditional 2SLS SCPC | Strong-identification inference with the published qualification | Numerical rank does not establish instrument strength; under endogeneity the finite conditional Gaussian argument is heuristic |
| Levels SPUR | Low-frequency Gaussian benchmark test and applicable paper asymptotics | Geometry, process assumptions, and a suitable low-frequency dimension |
| Residual SPUR | Conditional/residual benchmark test under its design-disturbance assumptions | Independence as in MW2024 section 4.5, or an applicable asymptotic relaxation; normal equations or E[u given X]=0 alone do not establish the required covariance law |
| Prespecified transformed regression | The paper's simulation-supported practical procedure | A general validity theorem for arbitrary integrated regressors and disturbances is not claimed |
| SPUR-selected transformation/pipeline | Published applied workflow with reported diagnostics | Selection adds an additional issue; no uniform post-selection coverage claim |
| Multiple responses or coefficients | Separate marginal inference | No familywise, simultaneous, or cointegrating-rank claim |

Sources: [MW2022](sources.md#mw2022), [MW2023](sources.md#mw2023),
[MW2024](sources.md#mw2024), and [BBV2026](sources.md#bbv2026). **[P, R, D]**

The Gaussian benchmark, the relevant asymptotic argument, the finite numerical
approximation, and empirical performance MUST remain distinguishable in
documentation and result summaries.

Positive numerical rank permits a calculation; it does not prove an asymptotic
approximation for q growing to the full available rank. Default SPUR q=15 is a
paper recommendation. The half-life q<=100 limit is a workload/support choice.
When q>=40 or all available positive directions are used, report an explicit
large-q diagnostic explaining the fixed-low-frequency limitation. **[R, D]**
Users MUST choose q before examining the tested response, or by a documented
geometry/design-only rule. The built-in SCPC selection rule meets that
restriction. Looking for a favorable p-value over q values is outside the
reported inference.

The clustered conditional model in scpc.md fixes a particular rank-one
observation pattern within each declared site. It MUST NOT be described as
protection against every possible within-cluster covariance. Correct multiway-FE
projection likewise does not establish asymptotic validity under arbitrary
growth of the number of effects or concentrated leverage. No arbitrary leverage
cutoff is imposed by this specification.

Ordinary 2SLS intervals require strong identification. Scalar Anderson–Rubin
testing of structural null hypotheses is the appropriate separate route for
weak-instrument robustness; support for general joint or weak-IV confidence
regions is not implied by supporting the 2SLS fit.

### 3.3 Constant span

SCPC accepts an identified supported regression with or without an intercept.
Residual SPUR requires that the full projected OLS design span the constant
vector. This makes the Brownian representation -M_X D M_X/2 valid after
projection by cancelling origin terms. **[I, D]** The requirement MUST be checked
on the actual design space, not merely by finding a column named "Intercept".
It is a package support restriction, not a claim that every alternative
no-intercept formulation is theoretically impossible.

The combined transformed-regression workflow additionally requires one explicit
nonzero constant column, of any name or scale, that can be removed separately
from the slope columns. A full-indicator/no-separate-constant parameterization
can span the constant but becomes rank-deficient after whitening. Reject that
combined fit with a message to refit using an intercept and explicit contrasts;
do not invent an unreported coefficient reparameterization. Residual tests
retain the broader constant-span contract, and standalone transformations may
return rank-reduced numeric columns.

## 4. Fitted samples, formulas, and coefficient identity

### 4.1 Fitted-state contract

An adapter MUST provide the actual coefficient vector and names, structural
residuals, fitted design, instruments where applicable, coefficient score
directions or equivalent sufficient matrices, and FE residualization operator.
Bread uses total cross products, without a residual-variance or degrees-of-freedom
factor. Definitions and dimensions are in scpc.md. **[D, I]**

Keep fitted observations in their fitted order. All sample selections, subsets,
split selections, omitted observations, and host-specific removals MUST be
applied before matching external coordinates or site labels. R fixest extraction
must include all applicable obsRemoved/split selections after subset. Stata
must verify structural residual sum of squares and appropriate regressor or
instrument normal equations against the active fit.

External data-frame matching MUST use identifying row information. A shared
constant, a non-unique set of regressor values, or matching row count is not
verification of observation identity. An API accepting an explicitly prealigned
array may treat it as the caller's declaration of fitted-row order; it MUST
distinguish that contract from automatically verified frame matching.
Ambiguous matches require an error naming the required row identifier or
alignment mechanism.

Do not silently discard fitted observations with missing coordinates or site
labels. Doing so changes the fitted estimator. Reject the call and identify the
missing fitted rows. An absorbed or aliased requested coefficient is
non-estimable and MUST NOT receive a numerical interval.

### 4.2 Formula and pipeline sample

Standalone formula transformations MUST use the host's ordinary formula
evaluator for supported terms and expressions, reject offsets, evaluate on
complete observations, and pad excluded output rows with NaN/NA or the host
equivalent. Formula labels are not literal data-column names.

The combined workflow MUST establish one analysis sample before levels
diagnostics, residual diagnostics, transformations, and adjusted regressions.
It must contain complete evaluated responses/regressors, coordinates and
required grouping variables. Stages MUST NOT silently use different samples.
Report the sample size, excluded-row count and row identity. If a stage needs
an additional exclusion, stop or restart the entire requested workflow with an
explicitly reported new sample.

Transform every fitted numeric column, not the formula syntax. Remove the
separable constant column from the adjusted fit and verify that the remaining
transformed slope design is identified. The LBM-GLS adjusted regression has no
intercept because centered whitening annihilates the constant. Do not report
an estimated adjusted intercept or restore one by an undocumented correction.
The workflow rules are in
[half-life and transformations](half-life-and-transforms.md).

## 5. Spatial sites and distance geometry

### 5.1 Sites and duplicates

Without a site label, A is the identity and m=n. With a site label, each
observation belongs to exactly one site and all observations in that site
MUST have identical canonical coordinates. Use coefficient-specific sums of
observation influences at the site; do not average observations or refit an
aggregated model. Site-level covariance calibration gives equal weight to
declared off-diagonal site pairs when defining average correlation. **[P, D]**

Distinct declared sites may have identical coordinates. They remain distinct
regression sites. Duplicate coordinates MUST NOT be automatically merged,
jittered, or treated as independent spatial locations. Define

$$
f=\frac{\#\{(i,j):i\ne j,\ d_{ij}=0\}}{m(m-1)},\qquad
K_{ij}(\infty)=\mathbf 1\{d_{ij}=0\}.
$$

The block limit follows from the exponential kernel. **[I]** A separate
idiosyncratic/nugget component is not included. Adding one would be a new
covariance-family decision requiring its own calibration.

Co-location means equality of canonical coordinates, with no fuzzy merging
radius. If representable distinct coordinates cannot yield a positive working
distance, report numerical geometry failure instead of silently reclassifying
them as duplicates. All-co-located data have zero diameter and provide no
spatial contrast; reject methods requiring that contrast.

### 5.2 Metrics and units

Euclidean distances use the supplied coordinate system and units. Translation
and a common positive rescaling may be used internally to avoid overflow.
Coordinatewise standardization is not a harmless common rescaling and MUST NOT
be applied silently.

Geographic coordinates use latitude and longitude in degrees, valid latitudes
in [-90,90], canonical equivalent longitudes, and one longitude at an exact
pole. The metric is the **great-circle arc** on a sphere of radius
6,371,008.8 meters. Compute degree differences before converting to angles,
and use a stable great-circle formula. Chordal distance is not an alternative
within this profile. **[A, D]**

| Operation | Working/reporting distance convention |
| --- | --- |
| Geographic SCPC | Arc distance divided by Earth circumference |
| Euclidean SCPC | Input units, or an explicitly recorded common working scale |
| SPUR covariance/basis | Distance divided by the sample diameter Delta |
| Normalized half-life | Fraction of the sample diameter |
| Unnormalized Euclidean half-life | Input coordinate units |
| Unnormalized geographic half-life | Meters |
| Geographic ISO radius | Meters |
| Euclidean ISO radius | Input coordinate units |

If d_work=d_input/s, then c_input=c_work/s and h_input=s*h_work.
Return the scale and units with any public decay parameter; a Nyström working
rate is not interchangeable with an unscaled input-coordinate rate.

Normalizing by diameter does not make arc and chordal geometries equivalent.
Spherical Gaussian benchmark calculations also do not automatically inherit
every Euclidean increasing-domain theorem.

## 6. Error, warning, and uninformative-result policy

Errors MUST identify the operation, offending quantity, applicable condition,
and a useful correction. Include counts or parameter values when available;
do not dump matrices or suggest arbitrary jitter.

| Situation | Required behavior |
| --- | --- |
| Unsupported estimator, weight, offset, FE representation, or ambiguous sample | Reject before numerical inference; name supported alternative/input contract |
| Unidentified/aliased requested coefficient | Error naming the coefficient |
| Unattainable average correlation | Error reporting requested value, f, and distinct-site count; suggest a valid bound or reconsideration of site definition |
| q exceeds positive rank or cuts a tied block | Error with requested q, available rank and eligible block boundaries |
| Zero/invalid response projection | Degeneracy error; a directional statistic is undefined |
| Zero SCPC estimated variance | Return explicit degenerate inference status and no ordinary t/p/CI claim, or fail that coefficient explicitly; do not manufacture significance |
| Valid nonzero response but no discrimination within the profile/tolerance | Conservative p=1, unbounded cutoffs/full half-life domain, and uninformative_within_profile status as defined by the method chapter |
| Materially indefinite covariance or failed factorization | Numerical/model error; no diagonal jitter or uninformative substitution |
| Unsatisfied integration/root accuracy at budget | Accuracy error or explicit incomplete diagnostic result without finalized inference |
| Incomplete half-life search | Retain evaluated point diagnostics with incomplete status; do not label the result a completed confidence set |
| Finite-mesh calibration / completed heuristic inversion | Return the specified approximation status and concise interpretation |
| Large q or unattained reference power target | Return the diagnostic with the achieved quantity and interpretive limit |

An uninformative full-domain result is deliberately conservative and does not
prove all continuous-parameter covariance shapes are identical. It MUST NOT
hide bad input, zero projected data, interrupted computation, or solver failure.
An interface MAY fail a multi-coefficient request atomically or return named
per-coefficient statuses, but it MUST never shift coefficient identities or
present failed entries as successful.

Warnings SHOULD be consolidated per result. Routine finite-mesh status belongs
in structured metadata and a concise result note, rather than repeated messages
for every evaluated covariance.

## 7. Result contract and reproducibility

Every successful result MUST identify the following information **where
applicable to the operation**:

1. Package/core version, implemented methods version, and numerical-profile ID.
2. Host, linear-algebra provider, precision, and effective thread setting.
3. Fitted/analysis sample, observation/site counts and coordinate units/scale.
4. Requested and effective q, validated rank/block eligibility information,
   exact or Nyström basis, and approximation diagnostics. A partial eigensolve
   may report a validated rank lower bound and an unknown total rank; metadata
   MUST NOT force an otherwise unnecessary full decomposition.
5. Actual confidence level or rejection sizes and the marginal interpretation.
6. Seed/RNG or supplied-draw provenance, draw counts and block allocation when
   simulations are used; landmark seed/sampler when applicable.
7. Calibration/search status, relevant meshes, convergence/error diagnostics,
   and uninformative or incomplete status.

These are logical fields. Host-native names and containers may differ if their
mapping is documented and no information is lost. The core mathematical result
must be available without parsing printed text. Bounds may be infinite; do not
replace true infinity with a finite display sentinel.

Default simulation seed is 42; default Nyström landmark seed is 1, with a
uint32 landmark seed domain. Simulation RNGs are local to the call and MUST NOT
alter the caller's global RNG state. A seed alone does not promise cross-language
or cross-version bitwise equality. Record the RNG/sampler identity and version.
Supplied, identical valid draws are the cross-host numerical comparison
mechanism. **[D]**

Distributional invariance under an orthogonal basis rotation differs from
fixed-draw reproducibility. The canonical orientation and complete-block rules
in numerics.md support row-permutation comparisons. Floating-point providers
need agree within stated tolerances, not bit for bit. Symmetries involving
indistinguishable duplicate rows must be assessed through projectors and
inference, without inventing persistent row identifiers.

## 8. Versioning and conformance

Adopted method definitions, benchmark families, simulation conventions, or
numerical profiles that can materially change results require a recorded
methods/profile version change. An equivalent solver, factorization, storage
layout or backend may satisfy the same profile if it meets the required
accuracy, diagnostics and conformance evidence.

Historical MATLAB/Stata output matching is a provenance check, not an override
of this specification. Do not add production compatibility modes for reversed
projections, finite infinity sentinels, old Brownian proxies, or the undeclared
Student floor. Preserve archived fixtures with explanations of intended changes.

Until the [migration gates](verification.md) pass, spec/VERSION, core
SPEC_VERSION and host compatibility checks MUST continue to identify the
implemented older method. They must move together when 0.6.0 is actually
implemented. Publishing this specification alone is not that event.
