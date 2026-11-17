# Methods specification 0.2.0

Version 0.2.0 inherits the mathematical definitions and independent reference calculations for [exact conditional OLS SCPC](../0.1.0/scpc-conditional-ols.md) and [SPUR half-life](../0.1.0/spur-half-life.md) from version 0.1.0 without changing those targets. It adds the required [IV projection](scpc-iv.md) and records the [adopted implementation conventions](implementation-conventions.md) for the clean shared-core packages.

The version 0.1.0 directory remains immutable as a record of the original audit, including the historical package defects and then-open product choices. Its historical descriptions do not define the corrected packages' defaults.

The IV extension follows Müller and Watson (2023), §5.3. Cluster and absorbed-effect mappings are included where needed to define the projection. This version does not give a complete theoretical specification for Nyström approximation, weighted fits, SPUR diagnostic tests, or spatial transformations.

Validation combines independent matrix identities, null rejection quadrature, the inherited 29 reference checks, and public native/fallback regression tests. A software parity test supports implementation consistency; it does not establish size control over a continuum of covariance parameters.
