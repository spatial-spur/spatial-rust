# Sources and provenance for 0.6.0

The [provenance tags](README.md#1-authority-and-document-map) distinguish
published definitions **[P]**, published recommendations **[R]**, author
software **[A]**, package decisions **[D]**, and algebraic consequences **[I]**.
This register identifies sources. It does not make historical numerical
conventions normative by reference.

## MW2022

Müller and Watson (2022), *Spatial Correlation Robust Inference*,
Econometrica 90(6), 2901–2935.

- [DOI 10.3982/ECTA19465](https://doi.org/10.3982/ECTA19465)
- [Author PDF](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf)
- Relevant locations: section 2, equations (3)–(6), Lemma 1, and Remarks 2.3
  and 2.6. These distinguish the covariance benchmark, Gaussian ratio
  calibration, expected-length objective, correlation-bound calibration, and
  recommendation to select q at 5%.

## MW2023

Müller and Watson (2023), *Spatial Correlation Robust Inference in Linear
Regression and Panel Models*, Journal of Business & Economic Statistics
41(4), 1050–1064.

- [Author PDF](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf)
- [Supplement](https://www.princeton.edu/~mwatson/papers/SpatialRegression_SupplementaryMaterial_20220907.pdf)
- Relevant locations: sections 2.1–2.3 and 5.1 for regression influence,
  sign/site benchmark and combined calibration; section 5.2/equation (10)
  for landmark approximation; section 5.3 for IV and its qualifications.
- The conditional weight-column order follows the actual residual operator.
  The transposed IV operator is algebraically required even when an older
  software implementation used another order.

## MW2024

Müller and Watson (2024), *Spatial Unit Roots and Spurious Regression*,
Econometrica 92(5), 1661–1695.

- [DOI 10.3982/ECTA21654](https://doi.org/10.3982/ECTA21654)
- [Author PDF](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf)
- Relevant locations: section 2.2/equation (2) for Lévy–Brownian covariance;
  section 2.5 for local-to-unity covariance; sections 4.1–4.4/equations
  (16)–(20) for low-frequency diagnostics and half-life inversion; section
  4.5/Remark 4.5 for residual-design assumptions; section 5/equation (21)
  for transformations and their simulation-supported use.
- A published tuning recommendation and a general validity theorem have
  different authority. The specification records that distinction.

## BBV2026

Becker, Boll and Voth (2026), *Testing and correcting for spatial unit roots in
regression analysis*, Stata Journal 26(2), 177–202.

- [DOI and publisher article](https://journals.sagepub.com/doi/10.1177/1536867X261449932)
- Relevant locations: section 2.7 for applied half-life interpretation,
  section 3 for transformations, section 4.1 for the practitioner decision
  tree, and Appendix A for empirical replication.
- This is provenance for an applied workflow, not a general uniform-coverage
  theorem for data-selected transformations.

## Author replication software

The SPUR Econometrica replication archive is
[Zenodo version 1, DOI 10.5281/zenodo.11199509](https://doi.org/10.5281/zenodo.11199509).
The inspected file is ReplicationPackage_MS21654_MuellerWatson_r2.zip.

SHA-256:
9c4348874ace003a855338653babc4202c2f8e8af32216b6f0c6af3a83bb8dd0

Paths below are relative to its ReplicationPackage_MS21654_MuellerWatson_r2
directory:

| File | Provenance relevance |
| --- | --- |
| cscpc/cscpc.m | Historical finite grid, q search, combined intervals and baseline p-value convention |
| cscpc/findcv.m | Student-t starting floor and upper-bracket numerical convention |
| cscpc/get_cv_eta.m | Conditional OLS projection order |
| cscpc/get_cv_eta_cluster.m | Site-normalized coefficient direction and full-design projection |
| cscpc/getcbar_W_large_n.m | Author landmark normalization/pooling profile |
| Matlab/functions/Spatial_I0_Test.m | I(0) reference experiment and normalization |
| Matlab/functions/Spatial_I1_Test.m | I(1) tuning experiment and normalization |
| Matlab/functions/getpow_qf.m | Simulated power calculation |
| Matlab/functions/c_ci.m | Historical Brownian proxy, half-life mixture and empirical inversion |
| Matlab/Chetty_Regressions_Table_1_ComputeResults.m | Replication draws, grids and hull |
| Matlab/Chetty_Regressions_Table_1_FormatTable.m | Finite-sentinel reporting convention |
| matlab_functions/lbm_gls_matrix.m | Inherited 1e-10 whitening cutoff at unit diameter |
| matlab_functions/getdistmat_lat_lon.m | Geographic distance implementation |

The September 2022 JBES replication separately contains Matlab/get_cv_eta.m
and Fortran/Sigmod.f90, corroborating the conditional projection order.
Its inspected Matlab/get_cv_eta.m SHA-256 is
99f001c89f8bb819685c4aa27abc69211fe6d273c5ee162b7533962b1ca8ae77.

The [historical source register](../0.1.0/sources.md) preserves additional
artifact digests and package pins. It is archival evidence, not a dependency
for the current method definition.

## Davies

Davies, author's implementation and documentation for linear combinations of
chi-squared random variables and related quadratic-form probabilities:

- [Author documentation and software](https://www.robertnz.net/QF.htm)
- [Algorithm paper](https://www.robertnz.net/pdf/lc_chisq.pdf)

This source motivates the optional deterministic SPUR solver investigation.
Accuracy controls and failure flags are part of such a solver's contract.
No performance benefit is asserted without measurement.

## Baseline and audit evidence

The reviewed implementation is
[spatial-rust at c41aaaca](https://github.com/spatial-spur/spatial-rust/tree/c41aaaca50559d2750bf10b2d5e8d2d1aac2746d).
It reports methods 0.5.0. Its code is implementation evidence, not authority
over this adopted 0.6.0 specification.

The 24 November 2026 methods audit established:

- A duplicate-site SCPC Student-floor/p-value discrepancy larger than the
  declared cutoff tolerance, conservative over the evaluated family.
- Accurate ten-node half-life mixture integration on the tested geometries,
  without claiming a universal integration bound.
- A costly but explicitly finite half-life refinement policy, without a newly
  demonstrated missed actual confidence component.
- Missing definitions, provenance labels, and numerical-scope qualifications.

The immediately preceding review of the same commit supplied the continuous
tail counterexample and timing/memory comparisons. They motivate verification
requirements; they are not proofs of general performance or coverage.
The [verification chapter](verification.md) defines the evidence an
implementation must supply to declare conformance.
