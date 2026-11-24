# Verification and migration to methods 0.6.0

This is the implementation acceptance plan for the adopted specification.
It distinguishes documentation completion from numerical implementation.
The reviewed baseline at c41aaaca reports methods 0.5.0. **None of the new
requirements below is deemed implemented merely because this chapter exists.**

## 1. Conformance rule

Before a binary advertises methods 0.6.0:

1. Implement every applicable normative requirement in this directory.
2. Supply the independent checks and host validations below, with reproducible
   fixtures, dependencies, settings and numerical-profile metadata.
3. Record intended differences from 0.5.0 and author results.
4. Update spec/VERSION, core SPEC_VERSION, host version guards and
   installed-artifact expectations together.

A missing mandatory host capability must be removed from the claimed support
matrix for that release or block the relevant host release. A successful Python
run is not evidence that R or Stata extraction and presentation conform.
No package version, backend choice or documentation heading substitutes for
the implemented-method identifier.

## 2. Baseline changes required

| Area | Required change relative to the reviewed baseline | Completion evidence |
| --- | --- | --- |
| SCPC Student floor | Permit downward critical-value bracketing; no undeclared IID addition | V01, V02 |
| SCPC finite tail | Add every prescribed finite-tail point, use common mesh for q/CV/p, expose actual rates/status | V03 |
| I(0) extreme tail | Replace silently coarsened wide/dense fallback with the declared budget error | V04 |
| Probability convergence | Separate depth/iteration exhaustion from satisfied accuracy; error-aware bracket checks | V02, V05 |
| Decay/root and numerical metadata | Report attained criterion/bracket; no silent looser fallback | V05 |
| Fitted-state support | Reject unsupported estimator subclasses and missing/misaligned fitted state | V10 |
| SPUR diagnostics | Canonical reference_tuning_power, precise fallback/status and normalization metadata | V06 |
| Uninformative cases | Conservative p=1/full-domain outputs for valid nonzero contrasts under the profile; no substitution for failure | V07 |
| Half-life output | Reusable point-null evaluation; all actual evaluations/counts; mapped transition brackets; completed heuristic versus incomplete status | V08 |
| Efficiency | Reuse geometry/decay limits through q expansion and geometry/factors/point evaluations within half-life | V11 |
| Combined workflow | Require separable constant column; preserve one sample and identifiable transformed slopes | V10 |
| Documentation/results | Scope, provenance, metric/units, numerical domains and applicable diagnostics consistent across hosts | V09, V10, V12 |

Some requirements restore or clarify behavior already present in parts of the
baseline. Verify that behavior rather than rewriting it unnecessarily.

## 3. Numerical and algebraic verification

### V01. Common-family p-values and intervals

Use a public fitted-model call with n=100, an intercept, 40 coordinates equal
to zero, and remaining coordinates 1,...,60. Use a deterministic response
generated from a recorded standard-normal seed. The co-location fraction is
40*39/(100*99), and include avc=f+0.001.

At the baseline, q=4, the inherited critical value was 2.7764451052 and the
independently required cutoff on its old evaluated grid was 2.7655437719.
A constant shift produced p=0.0496195862 while the reported 95% interval
contained zero. These numbers identify the original defect; they are **not**
the new profile's golden cutoff because its mesh is different.

For 0.6.0:

- Construct the actual full declared mesh and each relevant projected family
  independently from the returned basis and the mathematical formulas.
- Independently integrate fixed-covariance Gaussian probabilities.
- Check the returned critical bracket, the absence of an extra Student floor,
  and agreement of reported p and CI under the stated rounding/error allowance.
- Repeat at f+1e-8 and f+0.01; include default-avc controls with fewer
  duplicates, where no extra floor should have been active.

The test must compare the actual declared family, not a copied production
helper or an assumed IID endpoint.

### V02. Probability boundaries and conditional projections

Verify OLS FWL/bread scaling against a dense independent calculation. Check
that sign multiplication and residualization are in the required order.
For IV, use the actual transpose and verify reduction to OLS when P_Z X_w=X_w.
Include site aggregation and one-/multiway FE examples with observation-level
projection before aggregation.

Include singular PSD laws with a genuine finite atom, or the declared clear
unsupported-degeneracy error. Distinguish inclusive observed p-values from
strict rejection probabilities; do not assert equality of those functions
at an atom. Check zero variance, zero influence directions, no-intercept SCPC
and the prohibition on manufactured significance.

An error-aware critical bracket must have satisfied the same complete family
on both sides. Ordinary continuous-law p/CI differences must be explainable
by the recorded t bracket or assessed probability uncertainty.

### V03. Mesh, co-location and omitted-tail regression

Verify that exact and Nyström calls use full-site decay limits and that every
q candidate, cutoff, p-value and confidence level includes all prescribed
finite points plus K(infinity). Check the algebraic co-location endpoint
against summed-row outer products.

Retain the previous n=500/q=8 concentrated-leverage fixture as a diagnostic.
At the baseline, one omitted decay 44603.48665508749 had rejection probability
0.05098875005198544 at the returned cutoff 4.090183822704166. The new finite
tail must extend across that previously omitted scale range. Independently
verify its evaluated-family calibration and report the probability at that
off-mesh probe and additional nearby probes.

Report any remaining off-mesh difference explicitly; this verification does
not transform a finite mesh into a continuum certificate. Do not relabel
an off-mesh diagnostic as a formal uniform bound. Lower-leverage controls
help distinguish a targeted fixture from a claim about common default behavior.

### V04. Correlation-floor and tail-budget cases

Check avc domain endpoints and f<avc. For I(0), verify f<0.03, exact
f=0.03 limit-only inference, and f>0.03 infeasibility. Exercise the .001
reference-correlation fallback and the .0001 hybrid-grid target fallback.

Verify actual distinct mesh counts and endpoint inclusion, not a literal
"80" count. Check log-domain generation when a direct rate ratio overflows.
When more than 256 tail intervals or an unrepresentable tail endpoint is
required, assert the documented error. There must be no silent truncation,
coarsening, jitter or hidden fallback.

### V05. Accuracy and numerical failure

Compare Gaussian probabilities with an independent representation on realistic
well-scaled and moderately ill-conditioned covariance examples, including
co-location, mixed sign conditional weights and concentrated leverage. Record
absolute/relative errors in the probability or complement being tested.
Separate numerical integration error from a Monte Carlo discrepancy.

Force integration/root budget exhaustion in a controlled fixture and verify
accuracy failure rather than successful status. Test scale-aware PSD checks,
especially a materially indefinite small-variance block hidden beside a large
variance. Test stable coordinate/response scaling and valid small positive
distances without broadening the supported model.

Retain regressions for bugs already identified; do not add expensive sweeps
solely to enumerate unrealistic floating-point inputs.

### V06. SPUR definition, normalization and simulation ranks

Construct Brownian and residual projection matrices independently. Verify
I(0)/I(1) statistics, reference alternatives and the finite-rate raw
normalization maps. At the co-location fallback, verify its direct definition
rather than substituting an infinite number into a finite-rate formula.

Use fixed supplied independent draw blocks to compare hosts. Verify inclusive
plus-one p-values, strict cutoffs, exact integer order-statistic ranks and ties,
including B too small to reject at a requested size. Check reference-tuning
status and the meaning of achieved power; the final composite cutoff must
not be confused with the reference cutoff.

Check that user-supplied invalid norms/shapes fail and that the documented
distribution/independence requirement is visible. Batch results remain
marginal scalar tests; no joint or cointegrating-rank interpretation is added.

### V07. Uninformative versus undefined

Use a valid nonzero response and geometry-only covariance-shape equivalence
case. Verify conservative p=1, infinite thresholds or full half-life domain,
and uninformative_within_profile. No fake tuning success or simulation counts
may be fabricated.

Separately test zero projected response, invalid/indefinite covariance,
unsupported rank, and interrupted search. Those failures must not be
converted to an uninformative successful result.

### V08. Half-life mixture and inversion

Verify the angular density on a unit direction and its invariance to positive
covariance scaling. Independently integrate the uniform-in-h **density**
mixture, including determinant weights. Check true zero/infinite limits and
physical versus normalized units.

The audit compared 10-node with 40-node panels on a regular n=128/q=100 line,
an irregular n=80/q=15 geometry, four tight clusters, and a five-point close-pair
geometry. Its maximum relative density differences were approximately
1.10e-11, 4.88e-15, 7.11e-15 and 4.55e-13. Retain those independent
convergence checks or equivalent documented fixtures. They support the
retained order; they are not a uniform bound on every geometry or endpoint.

Check point-null acceptance p>alpha under the declared represented-level
boundary convention, all detected components, empty sets,
unbounded components and separately labelled hulls. Verify every actually
evaluated grid/refinement/bisection point and its count is available, with
phase labels and no invented counts for bypassed work. Check mapped endpoint
brackets and outward rounding.

Synthetic narrow gaps/islands may test search mechanics, but must not be cited
as a proof of actual-likelihood completeness. Check that all component results
carry heuristic scope and incomplete searches do not masquerade as completed
confidence sets.

### V09. Subspaces, geometry and transformations

Verify complete-block behavior under simple and tied eigenvalues, using
projectors for orientation-independent comparisons. Test row permutations,
translations and common unit changes within representable input precision.
For finite supplied simulations, distinguish canonical fixed-draw comparisons
from arbitrary rotated-basis distributional invariance.

Check full-landmark/exact subspace agreement and recorded seed/subset reuse
through q expansion. Verify pooled Frobenius weighting and post-Gram
orthogonality, without demanding a full eigensolve solely for diagnostics.

For transformations check Brownian pseudoinverse and constant annihilation,
the inherited strict whitening threshold, nearest positive distance and tied
neighbors, strict ISO boundaries/isolated rows, and cluster demeaning.
Geographic tests must distinguish arc from chordal distance and verify units.

### V10. Host/sample contract and Chetty comparison

Exercise each claimed host on fitted samples with exclusions, formula
expressions/interactions, duplicate row labels, FE selections, coefficient
aliases and ordinary IV where supported. Verify both correctly aligned arrays
and clear rejection of ambiguous frame matching or unsupported estimators.

The combined workflow must use one sample, actual evaluated design columns,
a separable constant, and an identified transformed slope fit with no added
intercept. A full-indicator design without a separable constant receives the
documented error; residual diagnostics can still accept its constant span.

Run the Chetty example across supported hosts with common supplied draws and
coordinate conventions. Compare point estimates, q, SEs, cutoffs, p-values,
half-life components and transformation output. Explain intended differences
from author/main fixtures due to corrected projections, exact Brownian
geometry, finite-simulation ranks, true endpoints and the adopted 0.6.0 mesh.
Agreement with an obsolete approximation is not a conformance requirement.

### V11. Efficiency checks

Instrument geometry/decay preparation and half-life point evaluation to verify
per-call reuse. A repeated q expansion must not recompute full-site decay roots;
repeated half-life points must not rerun their whole simulation.

Record representative exact/Nyström and half-life time and peak memory against
the prior implementation and main. Specify B explicitly: a 1,000-draw run
does not describe the 100,000-draw default. Track output differences alongside
speed. The auto threshold is evaluated as a memory/runtime policy, not assumed
to be a universal crossover.

No universal speedup or bitwise-provider requirement is imposed. Broaden
benchmarking only to resolve a concrete remaining performance risk.

### V12. Documentation and output consistency

Validate local links/anchors, source identifiers, equations, units, table
references and all D01–D20 mappings. Check that every delegated numerical
rule has a definition in the current chapter set.

Verify applicable result metadata, rank lower bounds where a partial solve is
used, stage/status names, implemented-method version, and clear warning/error
messages. Ensure no runtime output claims 0.6.0 before the implementation
gates pass.

## 4. Order of implementation

1. Implement the common SCPC family/bracketing and bounded tail policy, plus
   accuracy failure propagation. Verify V01–V05.
2. Add output/status contracts, conservative uninformative handling, complete
   half-life trace and point-null access; verify V06–V08.
3. Finish adapter/sample/estimator and transformation obligations and
   cross-host Chetty comparisons; verify V09–V10.
4. Verify per-call reuse and representative performance, resolve remaining
   mismatches, and update runtime method markers together; verify V11–V12.

This order is implementation guidance, not permission to ship a partially
conforming binary as methods 0.6.0. Documentation changes can be merged first
while the baseline continues to identify methods 0.5.0.
