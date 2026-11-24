# Adopted decisions and provenance

This register implements all 20 specification items from the 24 November 2026
audit. A decision labelled retained is still an explicit adopted choice.
The normative details are in the linked chapters; this document is the
traceability and rationale record.

## 1. Audit decision register

| ID | Adopted decision | Origin and rationale | Normative location |
| --- | --- | --- | --- |
| D01 | One complete versioned specification with an explicit incorporated chapter set; no implicit inheritance from historical versions | [D] Makes the current method independently implementable | [Authority](README.md#1-authority-and-document-map) |
| D02 | State separate Gaussian, asymptotic, empirical and numerical interpretations for each workflow | [P, R, D] A supported estimator or full rank does not establish scientific assumptions | [Scope](README.md#3-supported-calculations-and-scientific-assumptions) |
| D03 | Exact declared-site semantics; great-circle arc metric; explicit units/scales; true co-location endpoints; no automatic merging, jitter or nugget | [P, A, D, I] The chosen kernel determines its endpoint law | [Geometry](README.md#5-spatial-sites-and-distance-geometry) |
| D04 | Support actual unweighted OLS/ordinary 2SLS only as listed; reject penalized/weighted/offset and unsupported fits; record strict SCPC avc domain | [D] Clear unsupported-input errors are preferable to undeclared estimator substitutions | [Support](README.md#31-host-support), [decay domain](numerics.md#3-common-decay-calibration) |
| D05 | Restore coefficient bread, influence, SE scaling, OLS/IV transpose, FE projection and site aggregation formulas | [P, I] Required algebra, not an optional compatibility convention | [SCPC](scpc.md) |
| D06 | Remove Student-t lower floor; use one declared family for p and cutoffs; inclusive observed tails at atoms and strict cutoff/closed-CI boundary | [D, I] Prevents the reproduced duplicate-site inconsistency without adding an undeclared IID benchmark | [SCPC](scpc.md), [critical solver](numerics.md#6-gaussian-probabilities-and-scpc-critical-values) |
| D07 | Clear errors for invalid/undefined computations; explicit degenerate variance status; conservative uninformative outputs only for valid nonzero contrasts | [D] Avoids silent repairs and distinguishes absence of discrimination from numerical failure | [Error policy](README.md#6-error-warning-and-uninformative-result-policy) |
| D08 | Retain location-only 5% expected-length q objective and complete blocks; explicitly label bounded candidate search a heuristic | [P, R, A, D] The 5% rule is paper-recommended; global q optimality is not promised | [SCPC](scpc.md) |
| D09 | Retain inverse-root landmark extension, one Frobenius scale/sample and pooled leading subspace; full-site calibration; record sampler and memory purpose of auto threshold | [D] Declares the actual approximation rather than claiming exact paper-algorithm equivalence | [Nyström](scpc.md) |
| D10 | Adopt finite-mesh calibration explicitly; retain 1.05 SCPC interior, add common bounded 1.2 finite-tail probes plus exact limit; preserve I(0) hybrid interior | [D] Addresses the observed omitted tail without pretending to certify a continuum | [Meshes](numerics.md#4-scpc-covariance-mesh), [I(0)](numerics.md#5-i0-covariance-mesh) |
| D11 | Separate accuracy criteria and units; no silent success at depth/iteration/resource caps; remove silently coarsened tail fallback | [D] Termination is not convergence; rare unsupported numerical ranges get actionable errors | [Numerical profile](numerics.md) |
| D12 | Define raw SPUR statistic normalization and finite-rate paper/code parameter maps directly; define limit fallback separately | [A, D, I] Equivalent tests need not have identical raw statistics | [SPUR](spur.md) |
| D13 | Retain bounded reference tuning; use reference_tuning_power and explicit attainment/status; do not claim final composite power/global optimum | [R, A, D] The diagnostic reports the experiment actually performed | [SPUR](spur.md) |
| D14 | Independent Gaussian tuning/inference blocks; inclusive plus-one ranks; exact integer cutoffs; local RNG and supplied-draw assumptions | [D, I] Gives a precise finite-simulation contract and cross-host comparison mechanism | [SPUR](spur.md), [simulation profile](numerics.md#7-simulations-tuning-and-half-life-quadrature) |
| D15 | Exact Brownian basis and endpoint laws; angular-density mixture uniform in half-life; full null domain; components and separately labelled hull | [P, D, I] Preserves the intended mathematical target instead of replication proxies | [Half-life](half-life-and-transforms.md) |
| D16 | Retain bounded heuristic half-life search, expose point-null operation and every actual evaluation; report incomplete/heuristic status | [D] Completeness is not proved by more probes; no new certified optimizer is required | [Half-life](half-life-and-transforms.md) |
| D17 | Return transition brackets in compactified and half-life units; keep Monte Carlo resolution distinct from boundary digits | [D, I] Avoids an unsupported uniform precision claim in physical units | [Half-life](half-life-and-transforms.md) |
| D18 | Explicit LBM pseudoinverse and inherited cutoff; positive-distance NN ties; strict ISO radius and zero isolated rows; no adjusted intercept | [P, A, D] Resolves transformation semantics and their provenance | [Transformations](half-life-and-transforms.md) |
| D19 | One verified fitted/pipeline sample, actual formula evaluation, coefficient identity, conditional metadata, and qualified invariance guarantees | [D, I] Prevents host or row-order differences from changing the estimand silently | [Samples](README.md#4-fitted-samples-formulas-and-coefficient-identity), [results](README.md#7-result-contract-and-reproducibility) |
| D20 | Tag mathematical definitions, recommendations, author conventions and package choices; preserve archived fixtures without obsolete production modes | [D] Historical replication remains evidence, not a constraint on more faithful methods | [Sources](sources.md), [versioning](README.md#8-versioning-and-conformance) |

## 2. Important provenance distinctions

| Rule | What comes from a paper | What is adopted by this package |
| --- | --- | --- |
| SCPC inference | Exponential benchmark, spatial PC statistic and nuisance supremum | Finite mesh, full-site tail probes, solver tolerances and diagnostics |
| q objective | IID expected-length objective; recommendation to choose q at 5% | Initial candidate caps, boundary expansion, numerical blocks and stopping rule |
| Conditional inference | Actual residual operator and site/sign benchmark | Robust fitted-state extraction and explicit numerical degeneracy policy |
| Nyström | Landmark approximation and pooling concept | Inverse-root weighting, Frobenius scales, global pooled truncation, full-site calibration and auto threshold |
| I(0)/I(1) | Low-frequency benchmark tests | Draw blocks, ranks, finite mesh, raw normalization and attained-power reporting |
| I(0) tuning values | 0.03 bound, 0.001 reference, approximately 50% reference power are recommendations | Their exact finite computation, floor fallback and unattained-target handling |
| Half-life | Uniform-in-half-life angular-density mixture and test inversion | Ten-node dyadic integration, finite probes/refinement, complete diagnostics and search qualification |
| Brownian whitening | Symmetric inverse square root on the contrast space | Strict 1e-10 unit-diameter retention inherited from author code |
| ISO | Strict-radius neighborhood and zero isolated-row convention | Metric index and stable local difference accumulation |
| NN | Nearest-neighbor transformation idea | Nearest strictly positive distance and averaging numerical ties |
| Geographic application | Applicable Gaussian spatial benchmark | Arc metric, sphere radius, units, and explicit qualification of asymptotic scope |
| Practitioner pipeline | Published applied recommendation and simulations | Supported sample/formula contract and no general selection guarantee |

Complete numerical blocks and exact endpoint laws are adopted to implement
well-defined mathematical objects. They are not compatibility requirements
merely because they can change historical output.

## 3. Deliberate simplifications

1. **One inference family.** There is no separate undocumented IID floor,
   conditional-p-value shortcut, or legacy critical-value mode.
2. **One bounded tail generator.** SCPC and I(0) share a log-spaced finite tail
   and its error semantics. Extreme ranges fail clearly instead of activating
   a second wide/dense/coarsened algorithm.
3. **One current specification.** Readers need not resolve old version fragments.
4. **Per-call reuse.** Geometry, decay limits, factors and point evaluations are
   reused without requiring a persistent caching framework.
5. **Qualified search.** The package does not add a general certified optimizer
   to support a completeness claim it has not established.
6. **Applicability-aware diagnostics.** Transformations do not need fictitious
   q/CI fields, and partial eigensolves need not compute all eigenvalues solely
   to populate a rank label.
7. **Explicit intercept support.** The combined transformed fit requires a
   separable constant column and identifiable transformed slopes. General
   constant-span residual diagnostics remain supported; automatic coefficient
   contrast reparameterization is outside the combined workflow.
8. **Clear limits instead of repairs.** No coordinate jitter, silent ridge,
   outcome-selected q, hidden row dropping, or estimator substitution.

## 4. Choices retained pending measured evidence

The ten-node half-life mixture, bounded current search recipe, 5% q reference,
whitening cutoff, landmark defaults and default Monte Carlo inference remain
adopted. The audit did not justify replacing them wholesale. Their provenance,
scope and diagnostics are now explicit.

The general deterministic quadratic-form SPUR solver is an optional follow-up
experiment with a promotion gate, not an unimplemented mandatory branch of
0.6.0. Likewise, a continuous nuisance or complete-set certificate is outside
this profile. Future claims require new numerical evidence and an explicit
profile change.

## 5. Baseline migration

These decisions are normative targets, not claims that the reviewed 0.5.0
implementation already satisfies them. [Verification and migration](verification.md)
records required code/output changes and the evidence needed before the
implemented methods version advances. Documentation-only adoption does not
make old outputs 0.6.0 outputs.
