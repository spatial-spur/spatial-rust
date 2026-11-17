# Sources and provenance

Methods specification 0.1.0. Published definitions, original software and empirical replication, and package comparison pins have separate authority.

## Published methods

| Source | Relevant sections |
| --- | --- |
| Müller and Watson (2022), *Spatial Correlation Robust Inference*, Econometrica 90(6), 2901–2935. [DOI 10.3982/ECTA19465](https://doi.org/10.3982/ECTA19465) · [Author PDF](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf) | §2, equations (3)–(6), centered basis, statistic, size control, and i.i.d. expected-length criterion. Remark 2.3 calibrates the correlation bound. Lemma 1 gives a rejection-probability representation. |
| Müller and Watson (2023), *Spatial Correlation Robust Inference in Linear Regression and Panel Models*, JBES 41(4), 1050–1064. [Author PDF](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf) · [Supplement](https://www.princeton.edu/~mwatson/papers/SpatialRegression_SupplementaryMaterial_20220907.pdf) | §§2.1–2.3 and 5.1, OLS scores, conditional sign model, conditional and location cutoffs, and their covariance calculation. Equations (3), (5), (7), (8). The factor $k^2/q$ in this spec follows directly from the rejection inequality in equation (8). |
| Müller and Watson (2024), *Spatial Unit Roots and Spurious Regression*, Econometrica 92(5), 1661–1695. [DOI 10.3982/ECTA21654](https://doi.org/10.3982/ECTA21654) · [Author PDF](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf) | §2.2 equation (2) gives the Lévy–Brownian covariance. §2.5 defines the local-to-unity model. §4.1 sets up the low-frequency basis and §4.4 equation (20) gives the half-life score and inversion. Table 1 is an empirical comparison. |
| Becker, Boll and Voth (2026), *Testing and correcting for spatial unit roots in regression analysis*, Stata Journal 26(2), 177–202. [DOI 10.1177/1536867X261449932](https://doi.org/10.1177/1536867X261449932) | §2.7 gives the applied half-life interpretation and Stata usage. Appendix A compares the Stata results to the author empirical application. |

## Author software and replication

**JBES regression replication materials, September 2022.** The downloaded, extracted package was identified by `ReplicationFiles_SpatialRegression.txt` and checked against the [author's replication package listing](https://www.princeton.edu/~mwatson/publi). The following checksums identify the inspected source files.

| File | Reference | SHA-256 |
| --- | --- | --- |
| `Matlab/get_cv_eta.m` | Lines 8–15 multiply by the residualized regressor, project off the full design, and then multiply by sign. | `99f001c89f8bb819685c4aa27abc69211fe6d273c5ee162b7533962b1ca8ae77` |
| `Matlab/mc_unitline_model_3_dm.m` | Lines 62–74 use this function to obtain the conditional critical value for the paper's simulation. | `c0a5ef8707eaf51416190f58d5a54d63bb056af697a2887230f993688009bb38` |
| `Fortran/Sigmod.f90` | Lines 203–219 multiply the spatial weights by the regressor, residualize, and apply the cluster-normalized regressor. For a singleton with nonzero regressor this last factor is its sign. `Fortran/compute.f90:250–259` shows the singleton normalization. | `566456017a651abb2489fc056dcda4b17ffb8e9e80ed254ab05417c8b79a56da` |

**SPUR Econometrica replication archive, version 1.** [Zenodo DOI 10.5281/zenodo.11199509](https://doi.org/10.5281/zenodo.11199509). The file is `ReplicationPackage_MS21654_MuellerWatson_r2.zip`. Its verified SHA-256 is `9c4348874ace003a855338653babc4202c2f8e8af32216b6f0c6af3a83bb8dd0` and MD5 is `6929113ff1ff6f19d2a9fd1ff47341c8`. Paths below are relative to the top-level archive folder `ReplicationPackage_MS21654_MuellerWatson_r2/`.

| Member | Reference |
| --- | --- |
| `Matlab/Chetty_Regressions_Table_1_ComputeResults.m` | Lines 9, 20–32, 99–114 specify seed, draw count, grid, normalized distances, strict acceptance and hull. |
| `Matlab/functions/c_ci.m` | Lines 10–16 specify the exponential basis approximation. Lines 23–79 give the projected covariances, common draws recolored by the null, alternative density mixture, and empirical p-value. |
| `matlab_functions/getdistmat_normalized.m` | Lines 5–14 give the diameter normalization. |
| `Matlab/Chetty_Regressions_Table_1_FormatTable.m` | Lines 63–69 give the upper endpoint display convention. |
| `Matlab/mat/Table_1_Results.mat` and `Matlab/out/Table_1.tex` | Saved normalized half-life bounds and formatted Table 1. |
| `cscpc/get_cv_eta.m` | Lines 2–18 give the full projection and the order of the regressor and sign factors for unclustered C-SCPC. |
| `cscpc/get_cv_eta_cluster.m` | Lines 2–13 multiply by the cluster regressor matrix, residualize on the full design, and aggregate with the normalized cluster regressor matrix. |
| `cscpc/cscpc.m` | Lines 68–137 give the exact-mode grid and component search. Lines 145–166 give target-regressor partialling and combined critical values. Lines 180–209 give reported intervals and the baseline p-value. |
| `cscpc/form_omega_mat.m`, `cscpc/findqcv.m`, `cscpc/findcv.m` | Finite decay grid, component criterion and final bisection convention, respectively. |

**Author Stata SCPC.** [Pinned source `064a6227558cdbf9b206232378bd98e1f570524b`](https://github.com/ukmueller/SCPC/blob/064a6227558cdbf9b206232378bd98e1f570524b/src/scpc.ado). The inspected `src/scpc.ado` has SHA-256 `9a6eb3a4a206541dc0c1a1a1376339f1ac92b09c1e91875d68fa7f05ec8c6342` and matches the pinned local package source byte-for-byte. Lines 153–186 and 197–201 give bisection and component choice. Lines 575–580, 641–671, and 899–906 give the combined p-value and historical conditional operator order.

**Upstream Stata version history.** The public repository has one branch, `master`, with no tags or GitHub Releases. Its instructions install directly from `master/src`. Of the 30 public-history snapshots containing `src/scpc.ado` and `src/scpc.pkg` through `064a622`, 23 use unconditional inference only. The other seven contain the reversed order in the unclustered conditional path. Those seven snapshots have four distinct `scpc.ado` contents, at the commits linked below. The first adds conditional inference as the default and changes the help header to `1.1.0 Jan 2022`. The install manifest reports `1.0.0` throughout this history. [Upstream install instructions](https://github.com/ukmueller/SCPC#install) · [GitHub Releases](https://github.com/ukmueller/SCPC/releases) · [January 2022](https://github.com/ukmueller/SCPC/commit/9ff765b62b374392bbca71a06f31df8ea3e3bc2c) · [October 2022](https://github.com/ukmueller/SCPC/commit/7aa0abe27ab7e50855a29efcf410f07b205365f4) · [May 2025](https://github.com/ukmueller/SCPC/commit/45e48a13aeaaa81e015e0bc4d4df92a67a134d31) · [June 2025](https://github.com/ukmueller/SCPC/commit/064a6227558cdbf9b206232378bd98e1f570524b)

## Package comparison revisions

| Package profile | Pinned revision and source |
| --- | --- |
| Python SCPC release | `83716a88de329d834b7b981d9bab7d58791b1633`. [Conditional weights](https://github.com/spatial-spur/scpc-python/blob/83716a88de329d834b7b981d9bab7d58791b1633/src/scpc/utils/matrix.py#L341-L362) · [Inference](https://github.com/spatial-spur/scpc-python/blob/83716a88de329d834b7b981d9bab7d58791b1633/src/scpc/core.py#L293-L319) |
| R SCPC release | `e6e451cafacdfc6337e3ba398efc3eaecdf9b350`. [Conditional weights](https://github.com/spatial-spur/scpcR/blob/e6e451cafacdfc6337e3ba398efc3eaecdf9b350/R/utils-matrix.R#L166-L180) |
| Stata SCPC release | `064a6227558cdbf9b206232378bd98e1f570524b`. The author Stata source above is identical. |
| Milestone 1 Rust SCPC | `spatial-rust@123dfcb4b24b1155c1bf7e832cf3f285eca07524`, `crates/spatial-core/src/scpc_op.rs:298–345`. |
| Python SPUR release | `8c1a3b7b4f0151299f57d1283c2692cf92318ac5`. [Exponential basis and score](https://github.com/spatial-spur/spur-python/blob/8c1a3b7b4f0151299f57d1283c2692cf92318ac5/src/spur/utils/inference.py#L233-L317) · [Grid and inversion](https://github.com/spatial-spur/spur-python/blob/8c1a3b7b4f0151299f57d1283c2692cf92318ac5/src/spur/utils/inference.py#L321-L359) · [Public units and endpoint](https://github.com/spatial-spur/spur-python/blob/8c1a3b7b4f0151299f57d1283c2692cf92318ac5/src/spur/core.py#L874-L907) |
| R SPUR release | `ca70e08d1b80d2dc04205e654da6f7485e8af442`. [Persistence calculation](https://github.com/spatial-spur/spuR/blob/ca70e08d1b80d2dc04205e654da6f7485e8af442/R/spur_persistence.R) · [Public units and endpoint](https://github.com/spatial-spur/spuR/blob/ca70e08d1b80d2dc04205e654da6f7485e8af442/R/spur_formula_api.R#L549-L578) |
| Stata SPUR release | `faaac276080e6bf389705ac00742db5841bb93d8`. [Persistence calculation](https://github.com/spatial-spur/spur-stata/blob/faaac276080e6bf389705ac00742db5841bb93d8/mata/compute_spur_persistence.mata) · [Public units and endpoint](https://github.com/spatial-spur/spur-stata/blob/faaac276080e6bf389705ac00742db5841bb93d8/spurhalflife.ado#L27-L56) |
