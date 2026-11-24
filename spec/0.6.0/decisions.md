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
| D04 | Support actual unweighted OLS/ordinary 2SLS only as listed; reject penalized/weighted/offset and unsupported fits; use the mathematical SCPC domain f<avc<1 | [D] Actual precision/rank/resource checks replace arbitrary historical avc bounds | [Support](README.md#31-host-support), [decay domain](numerics.md#3-common-decay-calibration) |
| D05 | Restore coefficient bread, influence, SE scaling, OLS/IV transpose, FE projection and site aggregation formulas | [P, I] Required algebra, not an optional compatibility convention | [SCPC](scpc.md) |
| D06 | Remove Student-t lower floor; use one declared family for p and cutoffs; inclusive observed tails at atoms and strict cutoff/closed-CI boundary | [D, I] Prevents the reproduced duplicate-site inconsistency without adding an undeclared IID benchmark | [SCPC](scpc.md), [critical solver](numerics.md#6-gaussian-probabilities-and-scpc-critical-values) |
| D07 | Clear errors for invalid/undefined computations; explicit degenerate variance status; conservative uninformative outputs only for valid nonzero contrasts | [D] Avoids silent repairs and distinguishes absence of discrimination from numerical failure | [Error policy](README.md#6-error-warning-and-uninformative-result-policy) |
| D08 | Retain location-only 5% expected-length q objective and complete blocks; explicitly label bounded candidate search a heuristic | [P, R, A, D] The 5% rule is paper-recommended; global q optimality is not promised | [SCPC](scpc.md) |
| D09 | Retain inverse-root landmark extension, one Frobenius scale/sample and pooled leading subspace; full-site calibration; maintained unbiased landmark sampler instead of prescribed LCG | [D] Declares the actual approximation; the old bitstream has no mathematical authority | [Nyström](scpc.md), [sampler](numerics.md#23-scpc-candidate-expansion-and-landmark-sampler) |
| D10 | Adopt finite-mesh calibration explicitly; retain 1.05 SCPC interior, add common bounded 1.2 finite-tail probes plus exact limit; preserve I(0) hybrid interior | [D] Addresses the observed omitted tail without pretending to certify a continuum | [Meshes](numerics.md#4-scpc-covariance-mesh), [I(0)](numerics.md#5-i0-covariance-mesh) |
| D11 | Separate accuracy criteria and units; no silent success at depth/iteration/resource caps; remove silently coarsened tail fallback | [D] Termination is not convergence; rare unsupported numerical ranges get actionable errors | [Numerical profile](numerics.md) |
| D12 | Use unit-trace SPUR covariance shapes and a bounded alternative-strength parameter; retain paper/author normalization maps only as provenance | [D, I] Removes arbitrary covariance scales and unbounded multiplier bookkeeping without changing fixed-alternative test ordering | [SPUR](spur.md) |
| D13 | One bounded strength search with exact limits and shared rank-cutoff helper; report reference_tuning_power and explicit attainment/status | [R, D] Removes type-7 quantiles, separate 200-step searches, three correlation anchors and underflow-based stopping; no final composite power/global optimum claim | [SPUR](spur.md) |
| D14 | Independent Gaussian tuning/inference blocks with independent counts; inclusive plus-one ranks; exact integer cutoffs; explicit seeded local RNG or supplied draws | [D, I] Drops implicit global-RNG compatibility; aligned bases/coloring are required for identical finite-draw comparisons | [SPUR](spur.md), [simulation profile](numerics.md#7-simulations-tuning-and-half-life-quadrature) |
| D15 | Exact Brownian basis and endpoint laws; angular-density mixture uniform in half-life; full null domain; components and separately labelled hull | [P, D, I] Preserves the intended mathematical target instead of replication proxies | [Half-life](half-life-and-transforms.md) |
| D16 | Retain bounded heuristic half-life search, expose point-null operation and every actual evaluation; report incomplete/heuristic status | [D] Completeness is not proved by more probes; no new certified optimizer is required | [Half-life](half-life-and-transforms.md) |
| D17 | Return transition brackets in compactified and half-life units; separate exact endpoint membership from an infinite/zero outer bound | [D, I] An unresolved enclosure is not acceptance of the endpoint; Monte Carlo resolution differs from boundary digits | [Half-life](half-life-and-transforms.md) |
| D18 | Explicit LBM pseudoinverse and inherited cutoff; positive-distance NN ties; strict ISO radius and zero isolated rows; no adjusted intercept | [P, A, D] Resolves transformation semantics and their provenance | [Transformations](half-life-and-transforms.md) |
| D19 | One verified fitted/pipeline sample, actual formula evaluation, coefficient identity, conditional metadata, and projector-based invariance; no bespoke canonical eigenspace requirement | [D, I] Preserves the estimand and mathematical invariances without promising identical Monte Carlo ranks under arbitrary rotations | [Samples](README.md#4-fitted-samples-formulas-and-coefficient-identity), [results](README.md#7-result-contract-and-reproducibility) |
| D20 | Tag mathematical definitions, recommendations, author conventions and package choices; preserve archived fixtures without obsolete production modes | [D] Historical replication remains evidence, not a constraint on more faithful methods | [Sources](sources.md), [versioning](README.md#8-versioning-and-conformance) |

## 2. Important provenance distinctions

| Rule | What comes from a paper | What is adopted by this package |
| --- | --- | --- |
| SCPC inference | Exponential benchmark, spatial PC statistic and nuisance supremum | Finite mesh, full-site tail probes, solver tolerances and diagnostics |
| q objective | IID expected-length objective; recommendation to choose q at 5% | Initial candidate caps, boundary expansion, numerical blocks and stopping rule |
| Conditional inference | Actual residual operator and site/sign benchmark | Robust fitted-state extraction and explicit numerical degeneracy policy |
| Nyström | Landmark approximation and pooling concept | Inverse-root weighting, Frobenius scales, global pooled truncation, full-site calibration and auto threshold |
| I(0)/I(1) | Low-frequency benchmark tests | Draw blocks, ranks, finite mesh, unit-trace reporting, bounded alternative strength and attained-power reporting |
| I(0) tuning values | 0.03 bound, 0.001 reference, approximately 50% reference power are recommendations | Their exact finite computation, floor fallback and unattained-target handling |
| Half-life | Uniform-in-half-life angular-density mixture and test inversion | Ten-node dyadic integration, finite probes/refinement, complete diagnostics and search qualification |
| Brownian whitening | Symmetric inverse square root on the contrast space | Strict 1e-10 unit-diameter retention inherited from author code |
| ISO | Strict-radius neighborhood and zero isolated-row convention | Metric index and stable local difference accumulation |
| NN | Nearest-neighbor transformation idea | Nearest strictly positive distance and averaging ties under a distance-relative tolerance |
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
9. **One SPUR tuning calculation.** Shape normalization, a strength in [0,1],
   exact limiting covariances and the inference rank cutoff replace historical
   raw scales, type-7 quantiles and unbounded doubling/halving phases.
10. **Ordinary randomness and linear algebra.** A maintained RNG and validated
    eigenspaces replace a prescribed LCG and custom orientation machinery.
    Cross-host exact comparisons align actual subsets, bases and colored draws.
11. **Geometry without invented measurement uncertainty.** Nearest-neighbor
    ties depend on computed distances, not on the coordinate origin's ULP.
12. **One numerical implementation.** Hosts extract and present fitted state;
    the core owns the mathematics. Logical preparation/evaluation boundaries
    do not require a caching framework or a large public options surface.

## 4. Choices retained pending measured evidence

The ten-node half-life mixture, bounded current set-search recipe, 5% q
reference, whitening cutoff, landmark counts/auto threshold and default Monte
Carlo inference remain adopted. A simpler set search has not yet demonstrated
equivalent component detection; changing a whitening cutoff merely because it
is inherited is not an accuracy argument. Their provenance, scope and
diagnostics remain explicit. The reference-alternative tuning search and
landmark bitstream are deliberately revised as recorded above.

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

## 6. Follow-up review and design decisions

The follow-up review of the initial 0.6.0 draft identified the following gaps.
They are resolved in the normative chapters, not left as additional optional
compatibility modes.
The reviewed draft was
[d979f7d](https://github.com/spatial-spur/spatial-rust/commit/d979f7d00bea5f726bb126b8b7e1f5ee4aa4ed0a);
these revisions change the unimplemented target, not an already released
0.6.0 runtime.

| Review item | Assessment and resolution |
| --- | --- |
| Very short half-life work and intermediate decay overflow | Real numerical/resource gap, chiefly an extreme-range case. Bound mixture preparation to 64 panels and short-scale probes to 512; evaluate finite h directly, since overflow of log(2)/h need not invalidate exp(-log(2)d/h). |
| SCPC infinite cutoff | Real conditional-calibration degeneracy. Positive numerator with zero denominator variance implies cutoff infinity. Return a specific coefficient error before root finding; no retry loop or q reselection. |
| Infinite outer half-life bound | Real interpretation issue. Separate exact infinity acceptance from an unresolved outer enclosure, including a bracket ending at exactly rejected t=1. Apply the same principle to zero. |
| I(0) lowest log-grid offset | Wording defect: include delta_min explicitly in [delta_min,s), with s once. |
| R unseeded RNG behavior | Deliberate support decision: reject implicit seed=NULL and use explicit integer seeds or supplied draws. Do not preserve global-state consumption merely for compatibility. |
| I(1) underflow-based tuning stop | Removed with the bounded strength search and exact limiting covariance; no dependency on flush-to-zero behavior. |
| Singleton sign, local Nyström block scale and notation | Clarified locally and in the common notation contract. |

Two corrections to the review's four-site SCPC example matter for reproducible
verification. With coordinates (0,0,1,1), the off-diagonal co-location fraction
is **1/3**, not 1/6. With intercept and x=(1,3,-3,-1), its conditional denominator
weights sum to zero within both location groups, so the degeneracy holds at
**every decay**, not only infinity. The conclusion remains correct; V02 gives
the exact fitted-response fixture.

The extra nearest-neighbor defect is independently algebraic. For locations
(0,1,17/8), the point at 1 has distances 1 and 9/8. Adding exactly representable
10^15 to every coordinate leaves those distances unchanged but raises the old
coordinate-ULP allowance to about 1/4, admitting the farther neighbor. A
distance-relative 16 epsilon rule preserves the declared geometry. V09 retains
this regression case.

The revised SPUR normalization is an algebraic reparameterization at a fixed
alternative. Changing the bounded tuning search can select a different
alternative, and changing the sampler can select different landmark subsets;
those result changes are intentional. The implementation gates require
rank/ordering checks, reference-power comparisons and cross-host examples.
This document makes no unmeasured claim that the new tuning rule is faster or
more powerful on every geometry.
