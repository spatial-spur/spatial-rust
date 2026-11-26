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
| D08 | Retain location-only 5% expected-length q objective and complete blocks (under the unified 64 eps m_op rank/gap rule); simplify candidate search to initial prefix 20 and doubling 2k | [P, R, A, D] The 5% rule is paper-recommended; removing the 4-branch lookup preserves default behavior while simplifying growth; global q optimality is not promised | [SCPC](scpc.md), [blocks](numerics.md#21-rank-and-eigenvalue-blocks) |
| D09 | Retain inverse-root landmark extension, one Frobenius scale/sample and pooled leading subspace; full-site calibration; maintained unbiased landmark sampler instead of prescribed LCG | [D] Declares the actual approximation; the old bitstream has no mathematical authority | [Nyström](scpc.md), [sampler](numerics.md#23-scpc-candidate-expansion-and-landmark-sampler) |
| D10 | Adopt one single-solve uniform log-rate mesh for both SCPC and I(0): step ratio at most 1.10 from c0 to max(c0,36/d_min) (N<=512) plus exact K(infinity) | [D] Eliminates cmax, the second root solve, and the two-zone seam; unifies SCPC and I(0) on one simple grid gated by full I(0) validation checks | [Meshes](numerics.md#4-shared-covariance-mesh), [I(0)](numerics.md#5-i0-covariance-mesh) |
| D11 | Separate accuracy criteria and units; no silent success at depth/iteration/resource caps; remove silently coarsened tail fallback | [D] Termination is not convergence; rare unsupported numerical ranges get actionable errors | [Numerical profile](numerics.md) |
| D12 | Use unit-trace SPUR covariance shapes and a bounded alternative-strength parameter, anchoring I(1) decay by contrast scale c_ref=n/(2 tr(Omega_B)) as a tuning heuristic; retain paper/author normalization maps only as provenance | [D, I] Removes arbitrary covariance scales and unbounded multiplier bookkeeping while scaling multiscale/outlier-absorbed I(1) tuning without claiming universal interval interiority | [SPUR](spur.md) |
| D13 | One bounded strength search with exact limits and shared rank-cutoff helper; report reference_tuning_power and explicit attainment/status | [R, D] Removes type-7 quantiles, separate 200-step searches, three correlation anchors and underflow-based stopping; no final composite power/global optimum claim | [SPUR](spur.md) |
| D14 | Independent Gaussian tuning/inference blocks; symmetric Omega^{1/2} draw coloring; inclusive plus-one ranks; exact integer cutoffs; explicit seeded local RNG or supplied draws | [D, I] Unique symmetric square roots make one base-draw rotation E_new=O' E_old exact across all covariances simultaneously, giving exact joint finite-B rotation invariance | [SPUR](spur.md), [simulation profile](numerics.md#7-simulations-tuning-and-half-life-quadrature) |
| D15 | Exact Brownian basis and endpoint laws; angular-density mixture uniform in half-life; full null domain; components and separately labelled hull | [P, D, I] Preserves the intended mathematical target instead of replication proxies | [Half-life](half-life-and-transforms.md) |
| D16 | Single uniform log-h initial half-life grid (d_min/32 to 10^6 at 8 probes/doubling, exact endpoints 0 and infinity) plus 20-step bisection in t; expose point-null API and every evaluation | [D] Replaces stitched 3-part initial grid, deduplication, 3-round midpoint refinement, and tau_t with a single grid bounded automatically by the 64-panel preflight cap | [Half-life](half-life-and-transforms.md) |
| D17 | Return transition brackets in compactified and half-life units; separate exact endpoint membership from an infinite/zero outer bound | [D, I] An unresolved enclosure is not acceptance of the endpoint; Monte Carlo resolution differs from boundary digits | [Half-life](half-life-and-transforms.md) |
| D18 | Explicit LBM pseudoinverse and inherited cutoff; positive-distance NN ties; strict ISO radius and zero isolated rows; no adjusted intercept | [P, A, D] Resolves transformation semantics and their provenance | [Transformations](half-life-and-transforms.md) |
| D19 | One verified fitted/pipeline sample, actual formula evaluation, coefficient identity, conditional metadata, and projector/Omega^{1/2} invariance; no bespoke canonical eigenspace requirement | [D, I] Symmetric square-root coloring gives exact joint invariance across rotations without needing bespoke eigenspace anchors | [Samples](README.md#4-fitted-samples-formulas-and-coefficient-identity), [results](README.md#7-result-contract-and-reproducibility) |
| D20 | Tag mathematical definitions, recommendations, author conventions and package choices; preserve archived fixtures without obsolete production modes | [D] Historical replication remains evidence, not a constraint on more faithful methods | [Sources](sources.md), [versioning](README.md#8-versioning-and-conformance) |

## 2. Important provenance distinctions

| Rule | What comes from a paper | What is adopted by this package |
| --- | --- | --- |
| SCPC inference | Exponential benchmark, spatial PC statistic and nuisance supremum | Single uniform log-rate mesh, solver tolerances and diagnostics |
| q objective | IID expected-length objective; recommendation to choose q at 5% | Initial prefix 20, doubling 2k, unified 64 eps m_op blocks and stopping rule |
| Conditional inference | Actual residual operator and site/sign benchmark | Robust fitted-state extraction and explicit numerical degeneracy policy |
| Nyström | Landmark approximation and pooling concept | Inverse-root weighting, Frobenius scales, global pooled truncation, full-site calibration and auto threshold |
| I(0)/I(1) | Low-frequency benchmark tests | Draw blocks, Omega^{1/2} coloring, ranks, finite mesh, unit-trace reporting, bounded alternative strength and attained-power reporting |
| I(0) tuning values | 0.03 bound, 0.001 reference, approximately 50% reference power are recommendations | Their exact finite computation, floor fallback and unattained-target handling |
| Half-life | Uniform-in-half-life angular-density mixture and test inversion | Ten-node dyadic integration, single log-h grid, 20-step bisection, complete diagnostics and search qualification |
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
2. **One single-solve uniform log-rate mesh generator.** SCPC and I(0) share
   one one-solve log-rate mesh (`ratio <= 1.10`, `N <= 512`) up to
   `max(c0, 36/d_min)` plus exact `K(infinity)` and its preflight error
   semantics. There is no `cmax`, two-solve stitched seam, ~79-solve hybrid
   I(0) grid, or coarsened fallback.
3. **One unified spatial rank/block-gap rule.** Every $m_{\mathrm{op}}\times m_{\mathrm{op}}$
   spatial operator uses $\tau_{\mathrm{floor}}=64\epsilon m_{\mathrm{op}}$, positive
   threshold $\max(\tau_{\mathrm{floor}},10^{-12}\lambda_{\max})$, and block-gap
   threshold $\max(\tau_{\mathrm{floor}},10^{-11}\lambda_{\max})$, instead of
   five separate operator-specific thresholds.
4. **Simple SCPC prefix growth.** Start at `k = min(limit, 20)` and double
   `2k` instead of maintaining a 4-branch `avc` lookup table and `ceil(k/2)`.
5. **Exact joint rotation invariance via symmetric square roots.** Coloring
   simulation draws with $\Omega^{1/2}$ replaces lower-triangular Cholesky
   factors so a single draw rotation $E_{\mathrm{new}}=O^\top E_{\mathrm{old}}$
   couples every covariance simultaneously across I(1), I(0), tuning, and
   half-life search.
6. **Single uniform log-$h$ half-life grid and 20-step bisection.** One
   8-probes-per-doubling grid from `d_min/32` to `10^6` (plus exact `0` and
   `infinity`, bounded at `<=667` initial points by the 64-panel quadrature
   cap) and 20-step bisection in `t` replace the 3-part stitched grid,
   deduplication threshold, 3 midpoint rounds, and `tau_t`.
7. **One current specification.** Readers need not resolve old version fragments.
8. **Per-call reuse.** Geometry, decay limits, factors and point evaluations are
   reused without requiring a persistent caching framework.
9. **Qualified search.** The package does not add a general certified optimizer
   to support a completeness claim it has not established.
10. **Applicability-aware diagnostics.** Transformations do not need fictitious
    q/CI fields, and partial eigensolves need not compute all eigenvalues solely
    to populate a rank label.
11. **Explicit intercept support.** The combined transformed fit requires a
    separable constant column and identifiable transformed slopes. General
    constant-span residual diagnostics remain supported; automatic coefficient
    contrast reparameterization is outside the combined workflow.
12. **Clear limits instead of repairs.** No coordinate jitter, silent ridge,
    outcome-selected q, hidden row dropping, or estimator substitution.
13. **One SPUR tuning calculation.** Shape normalization, a strength in [0,1],
    exact limiting covariances and the inference rank cutoff replace historical
    raw scales, type-7 quantiles and unbounded doubling/halving phases.
14. **Ordinary randomness and linear algebra.** A maintained RNG, validated
    eigenspaces, and symmetric $\Omega^{1/2}$ coloring replace a prescribed LCG
    and bespoke orientation machinery. Cross-host exact comparisons align
    actual subsets and retained bases (or apply $E_{\mathrm{new}}=O^\top E_{\mathrm{old}}$).
15. **Geometry without invented measurement uncertainty.** Nearest-neighbor
    ties depend on computed distances, not on the coordinate origin's ULP.
16. **One numerical implementation.** Hosts extract and present fitted state;
    the core owns the mathematics. Logical preparation/evaluation boundaries
    do not require a caching framework or a large public options surface.

## 4. Choices retained pending measured evidence

The ten-node half-life mixture, 5% q reference, whitening cutoff, landmark
counts/auto threshold and default Monte Carlo inference remain adopted. Changing
a whitening cutoff merely because it is inherited is not an accuracy argument.
Their provenance, scope and diagnostics remain explicit. The single SCPC/I(0)
`1.10` log-grid (`N<=512`) and uniform log-$h$ half-life grid are gated on
independent validation checks in [verification](verification.md), while the
reference-alternative tuning search and landmark bitstream are deliberately
revised as recorded above.

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
| Very short half-life work and intermediate decay overflow | Real numerical/resource gap, chiefly an extreme-range case. Use one 64-panel (640-node) preflight cap (which bounds $h_{\min}=d_{\min}/32\ge2^{-63}$ and limits the single log-$h$ initial grid to <=667 points); evaluate finite $h$ directly, since overflow of $\log(2)/h$ need not invalidate $\exp(-\log(2)d/h)$. |
| SCPC infinite cutoff | Real conditional-calibration degeneracy. Positive numerator with zero denominator variance implies cutoff infinity. Return a specific coefficient error before root finding; no retry loop or q reselection. |
| Infinite outer half-life bound | Real interpretation issue. Separate exact infinity acceptance from an unresolved outer enclosure, including a bracket ending at exactly rejected t=1. Apply the same principle to zero. |
| Incomplete half-life search contract | Simplify host implementation: fail with an actionable error by default on interruption or budget/accuracy failure; optional `incomplete` partial diagnostics remain allowed (`MAY`) without resumable search state. |
| SCPC `avc` domain alignment & shared I(0) mesh | Align `scpc.md` with `numerics.md` (`f<avc<1`), and unify SCPC and I(0) onto one single-solve `1.10` log-rate mesh up to `max(c0,36/d_min)` (`N<=512`) plus `K(infinity)` instead of maintaining `cmax` or a separate ~79-solve hybrid grid. |
| I(1) multiscale/outlier compact decay | Anchor $c(w)=c_{\mathrm{ref}}\,w/(1-w)$ with $c_{\mathrm{ref}}=n/(2\operatorname{tr}(\Omega_B))$ as a geometry/contrast-adapted tuning heuristic so residual outlier-absorbed and multiscale geometries scale naturally without claiming universal interiority. |
| R unseeded RNG behavior | Deliberate support decision: reject implicit seed=NULL and use explicit integer seeds or supplied draws. Do not preserve global-state consumption merely for compatibility. |
| I(1) underflow-based tuning stop | Removed with the bounded strength search and exact limiting covariance; no dependency on flush-to-zero behavior. |
| Singleton sign, local Nyström rule, $\Omega^{1/2}$ rotation coupling, and notation | Unify spatial rank/block thresholds via $64\epsilon m_{\mathrm{op}}$; replace Cholesky draw coloring with unique symmetric $\Omega^{1/2}$ so single draw rotation $E_{\mathrm{new}}=O^\top E_{\mathrm{old}}$ gives exact joint finite-$B$ invariance across all covariances. |

Two corrections to the review's four-site SCPC example matter for reproducible
verification. With coordinates (0,0,1,1), the off-diagonal co-location fraction
is **1/3**, not 1/6. With intercept and x=(1,3,-3,-1), its conditional denominator
weights sum to zero within both location groups, so the degeneracy holds at
**every decay**, not only infinity. The conclusion remains correct; V02 gives
the exact fitted-response fixture.

The extra nearest-neighbor defect is independently algebraic. For locations
(0,1,17/8), the point at 1 has distances 1 and 9/8. Adding exactly representable
10^15 to every coordinate leaves those distances unchanged but raises the old
coordinate-ULP allowance to about 1/4, admitting the farther neighbor. A pure
distance-relative `1e-11` rule eliminates coordinate-origin ULP dependence while
still preserving intended ties on standard 0.1-degree (~1.4e-13) and 0.01-degree
(~1.4e-12) decimal geographic grids (whereas a 0.001-degree grid reaches
~1.4e-11 and exceeds `1e-11`). Large-origin rotated Euclidean grids whose
relative distance discrepancy exceeds `1e-11` are intentionally not
reconstructed. V09 retains both the translation and decimal-grid regression cases.

The revised SPUR normalization is an algebraic reparameterization at a fixed
alternative. Changing the bounded tuning search can select a different
alternative, and changing the sampler can select different landmark subsets;
those result changes are intentional. The implementation gates require
rank/ordering checks, reference-power comparisons and cross-host examples.
This document makes no unmeasured claim that the new tuning rule is faster or
more powerful on every geometry.
