# SCPC and SPUR methods specification

Current version **0.3.0**. See the [version 0.3.0 record](0.3.0/README.md) for sampled covariance validation, score signs, and the SPUR pipeline. Earlier records remain available.

The inherited theoretical scope covers exact conditional OLS SCPC for unclustered observations and half-life inference for a scalar spatial series. This specification defines the statistical targets and records the finite-sample conventions used to reproduce particular sources. It can be read and checked independently of Rust, the host packages, and their release fixtures.

| Document | Contents |
| --- | --- |
| [Current conventions](0.3.0/README.md) | Sampled covariance validation, unclustered score signs, SPUR pipeline samples and finite simulation. |
| [IV projection and version 0.2.0 conventions](0.2.0/README.md) | Required transposed IV weights and adopted validation, simulation, and unit policies. |
| [Conditional OLS SCPC](0.1.0/scpc-conditional-ols.md) | OLS coefficient and scores, spatial basis, recommended component selection, and control over both Gaussian benchmark families. |
| [SPUR half-life](0.1.0/spur-half-life.md) | Lévy–Brownian basis, angular Gaussian likelihood, uniform half-life mixture, and confidence-set inversion. |
| [Conventions and decisions](0.1.0/conventions.md) | Author replication and package behavior, and open choices for a reproducible product profile. |
| [Source register](0.1.0/sources.md) | Papers and equation locations, original source files, artifact digests, and package revisions. |
| [Reference calculations](0.1.0/reference/README.md) | Independent derivations, executable checks, and versioned expected values. |

## Use and versioning

- Apply the [conditional OLS projection order](0.1.0/scpc-conditional-ols.md#required-projection-order) in every package and backend, including the shared Rust core and host-language fallbacks.
- Cite the version and method section in implementation changes and tests. Compare a mathematical identity against the reference calculations and a release result against the corresponding compatibility fixture.
- Give result-changing updates to a method definition or adopted finite convention a new spec version. Retain prior version directories and record the source, rationale, affected fixtures, and method review with the change. A new solver, factorization, or storage layout can satisfy the same specification.
- Resolve choices marked open in [Conventions and decisions](0.1.0/conventions.md) before changing production defaults. Package parity remains tied to its named source revision. A disagreement with the published method belongs in the convention record with a reproducer and a decision.

Version 0.1.0 records the paper definitions and audits the profiles that inform implementation. It leaves product defaults open where author and package profiles differ. Clustered, weighted, absorbed-FE, and IV SCPC, Nyström SCPC, SPUR diagnostic tests, and transformations are outside this version.
