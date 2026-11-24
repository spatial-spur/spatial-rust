# SCPC and SPUR methods specification

[Version 0.6.0](0.6.0/README.md) is the adopted complete specification and next
implementation target. Its chapter map defines all normative material without
implicit inheritance from historical editions.

[Version 0.5.0](0.5.0/README.md) remains the implemented baseline identified by
spec/VERSION, the core and host compatibility guards. Those markers must advance
together only after the [0.6.0 conformance gates](0.6.0/verification.md) pass.
Adding the new specification alone does not claim a 0.6.0 implementation.
Earlier editions preserve historical algorithms and verification references.

| Document | Contents |
| --- | --- |
| [Version 0.6.0](0.6.0/README.md) | Complete normative scope and formulas; paper/replication/package provenance; common SCPC family without a Student floor; bounded finite tails; explicit numerical/search guarantees, diagnostics, and conformance gates. |
| [Version 0.5.0](0.5.0/README.md) | Finer SCPC/I(0) grids, 10-point GL & adaptive half-life search, full-site Nyström decay, aligned formulas, no adjusted intercept, R C-API & Stata normal-equation guards. |
| [Version 0.4.0](0.4.0/README.md) | Unified fitted-state, geometry, calibration, results, and execution contracts. |
| [Version 0.3.9](0.3.9/README.md) | ISO neighborhoods, centered LBM-GLS, conditional IV, and SPUR dimensions and tuning. |
| [Version 0.3.8 conventions](0.3.8/README.md) | Multiscale half-life inversion, response projection, attainable SPUR calibration, and fitted SCPC designs. |
| [Version 0.3.7 conventions](0.3.7/README.md) | Continuous half-life inversion and SPUR rank inference. |
| [Version 0.3.6 conventions](0.3.6/README.md) | SPUR residual response projection and native methods compatibility. |
| [Version 0.3.5 conventions](0.3.5/README.md) | Residual SPUR scope, half-life reporting, decay bracketing, SCPC variance and IV. |
| [Version 0.3.4 conventions](0.3.4/README.md) | SPUR neighbor and geographic conventions, intercept handling and SCPC model API. |
| [Version 0.3.3 conventions](0.3.3/README.md) | Gaussian tails, coordinate rescaling, varying slopes and component search. |
| [Version 0.3.2 conventions](0.3.2/README.md) | Design preparation, fixed-effect spacing, Gaussian boundaries and tied diagnostics. |
| [Version 0.3.1 conventions](0.3.1/README.md) | IV numerical rank and fixed-effect projection. |
| [Version 0.3.0 conventions](0.3.0/README.md) | Sampled covariance validation, unclustered score signs, SPUR pipeline samples and finite simulation. |
| [IV projection and version 0.2.0 conventions](0.2.0/README.md) | Required transposed IV weights and adopted validation, simulation, and unit policies. |
| [Conditional OLS SCPC](0.1.0/scpc-conditional-ols.md) | OLS coefficient and scores, spatial basis, recommended component selection, and control over both Gaussian benchmark families. |
| [SPUR half-life](0.1.0/spur-half-life.md) | Lévy–Brownian basis, angular Gaussian likelihood, uniform half-life mixture, and confidence-set inversion. |
| [Conventions and decisions](0.1.0/conventions.md) | Author replication and package behavior, and open choices for a reproducible product profile. |
| [Source register](0.1.0/sources.md) | Papers and equation locations, original source files, artifact digests, and package revisions. |
| [Reference calculations](0.1.0/reference/README.md) | Independent derivations, executable checks, and versioned expected values. |

## Use and versioning

- Implementing 0.6.0 requires the complete incorporated chapter set, including its [SCPC projection definitions](0.6.0/scpc.md) and [numerical profile](0.6.0/numerics.md).
- Cite the version and method section in implementation changes and tests. Compare a mathematical identity against the reference calculations and a release result against the corresponding compatibility fixture.
- Record result-changing method or numerical-profile updates with source, rationale, affected fixtures and verification. An equivalent solver, factorization or storage layout can satisfy the same specification.
- The [0.6.0 decision register](0.6.0/decisions.md) resolves its adopted defaults and supersedes earlier open-choice lists for that target. Historical output equality does not override corrected mathematical definitions.
- Keep implemented method markers honest while a specification is ahead of the code. The [migration plan](0.6.0/verification.md) records required changes and distinguishes them from documentation completion.

Version 0.1.0 records the paper definitions and audits the profiles that inform implementation. It leaves product defaults open where author and package profiles differ. Clustered, weighted, absorbed-FE, and IV SCPC, Nyström SCPC, SPUR diagnostic tests, and transformations are outside this version.
