# Milestone 1 Review Handoff: Shared Rust Core for Exact OLS SCPC

## 1. Repository Layout, Branches, and Commit SHAs

All legacy `dev-rust` work in the six upstream repositories was preserved intact via `git stash push -u -m "preserve dev-rust shared-core prototype"` before creating clean `rust-core` branches from the pinned release commits.

| Repository | Branch | HEAD Commit SHA | Base / Notes |
| --- | --- | --- | --- |
| `spatial-rust` | `main` | `123dfcb4b24b1155c1bf7e832cf3f285eca07524` | New Cargo workspace + pinned M1 fixtures (30 files) |
| `scpc-python` | `rust-core` | `9c9e0623d6954d8a55cb20a4e422d2d5ecc7b7b8` | Parent: `83716a88de329d834b7b981d9bab7d58791b1633` (`main`) |
| `scpcR` | `rust-core` | `cb35c7398361f817871745b6cbe067e46f818757` | Parent: `e6e451cafacdfc6337e3ba398efc3eaecdf9b350` (`main`) |
| `scpc-stata` | `rust-core` | `c1108e762f80e3e4206a77fdc45edfc21422d672` | Parent: `9a84d463780159449a869fa4db3c71616d8080aa` (strict-Mata fix); grandparent: `064a6227558cdbf9b206232378bd98e1f570524b` (`master`) |
| `spur-python` | `rust-core` | `8c1a3b7b4f0151299f57d1283c2692cf92318ac5` | Untouched in M1 (pinned at `main`) |
| `spuR` | `rust-core` | `ca70e08d1b80d2dc04205e654da6f7485e8af442` | Untouched in M1 (pinned at `main`) |
| `spur-stata` | `rust-core` | `faaac276080e6bf389705ac00742db5841bb93d8` | Untouched in M1 (pinned at `main`) |

Standalone strict-Mata prerequisite commit on `scpc-stata`:
- `9a84d463780159449a869fa4db3c71616d8080aa` (`fix(mata): declare local variables in setOmsWfin for strict Mata compatibility`), separating the baseline fix from the Rust integration commit.

## 2. Architecture and Whole-Operation Boundary

### Shared Numerical Core (`spatial-rust/crates/spatial-core`)
- Public Milestone 1 entry point: `spatial_core::run_exact_scpc(&ExactScpcRequest) -> Result<ExactScpcOutput, CoreError>` in `crates/spatial-core/src/scpc_op.rs` (re-exported in `crates/spatial-core/src/lib.rs`).
- Whole-operation design: hosts extract regression/design objects, enforce row alignment, and invoke Rust once for the numerical SCPC pipeline:
  - `ExactScpcRequest`: `coords` (`n x d`, row-major), `metric` (`Euclidean` or `GreatCircleFractionLonLat` / `GreatCircleFractionLatLon`), `avc`, `uncond`, `cvs`, `coef` (`k_use`), influence matrix `wj` (`n x k_use`), conditional direction `xj` and `xjs` (`n x k_use`), low-rank projection factors `proj_u` and `proj_v` (`n x proj_rank`, applying `A -> A - U (V' A)` internally), `cv_tolerance` (`0.001`), and `nthreads`.
  - `ExactScpcOutput`: `scpcstats` (`k_use x 6`: `[coef, se, tstat, pvalue, ci_lower, ci_upper]` with interval $\hat\beta \pm cv_{\max} SE$), optional `scpccvs` (`k_use x 4` at `32%`, `10%`, `5%`, `1%`), `wfin` (`n x (q + 1)`), `distmat` (`n x n`), `c0`, `cmax`, `cvfin`, `q`, and `tied`.
- Linear algebra backend (`crates/spatial-core/src/linalg.rs`):
  - Default Apple build links `Accelerate` (`cblas_dgemm`, `vvexp`, `dgeev_`).
  - `--features portable` and non-Apple targets use `matrixmultiply::dgemm` and a pure-Rust eigensolver.
- Thread control (`crates/spatial-core/src/threads.rs`):
  - `run_with_threads(Some(1), ...)` executes sequentially inside a cached single-thread Rayon pool without touching the global pool and **never** routes to interpreted fallback code. `backend` and `nthreads` are strictly orthogonal across all three hosts.

### Host Adapters and Eligibility in Milestone 1

Milestone 1 accelerates **exact** (`method != "approx"`, `neff < 4500` when `method == "auto"`), **unclustered**, **non-FE**, **non-IV** OLS fits. In `backend = "auto"` (and when tied eigenspaces occur with `tied == true` in `auto`), ineligible models route cleanly to the host fallback; in `backend = "rust"`, ineligible models error cleanly.

1. **Python (`spatial-rust/crates/spatial-py` + `scpc-python`)**
   - PyO3 module `spatial_spur_rust` (`API_VERSION == 1`) exposes `run_exact_scpc(...)`, releases the GIL via `py.allow_threads(...)`, validates finiteness/shapes/`avc`/`nthreads`, and normalizes row/column layout.
   - `scpc-python/src/scpc/types.py` adds `backend: str = "fallback"` to `SCPCResult`.
   - `scpc-python/src/scpc/core.py` adds `backend: str = "auto"` (`"auto" | "rust" | "fallback"`, honoring `SPATIAL_BACKEND`) and `nthreads: int | None = None`, plus optional dependency `rust = ["spatial-spur-rust>=0.2.0"]` in `pyproject.toml`.

2. **R (`spatial-rust/crates/spatial-r` + `scpcR`)**
   - `scpcR/src/init.c` registers pure-C `.Call` entry point `wrap__rust_run_exact_scpc` (13 `SEXP` args) with zero C++ dependencies (`Rcpp`/`RcppArmadillo`/`RcppParallel` are not used; `otool -L` confirms no `libc++` or `libstdc++`).
   - `scpcR/src/rust/` compiles static library `libscpc_rust.a` (referencing local workspace crates `../../../spatial-rust/crates/spatial-r` and `spatial-core`).
   - `scpcR/src/Makevars` and `scpcR/src/Makevars.win` automatically fall back to pure-R (`norust`) mode if `SPATIAL_RUST=0` or `cargo` is unavailable, tracking `.spatial_rust_mode` to invalidate `init.o` and `$(SHLIB)` across mode switches.
   - `scpcR/R/core.R` adds `backend = c("auto", "rust", "fallback")` (also honoring `getOption("spatial.backend")` and `Sys.getenv("SPATIAL_BACKEND")`) and `nthreads = NULL` (honoring `scpcR.threads`, `SCPCR_NUM_THREADS`, `mc.cores`, and `_R_CHECK_LIMIT_CORES_`).

3. **Stata (`spatial-rust/crates/spatial-stata` + `scpc-stata`)**
   - Builds `scpc-stata/src/scpc_rust.plugin` (tracked in `scpc-stata/src/scpc.pkg` as `f scpc_rust.plugin`; `scpc.pkg` retains `d Requires: Stata version 11`, while SPI v3 requires Stata 14.1+ and older Stata versions fall back cleanly to Mata if plugin loading fails).
   - Uses SPI v3 C bridge (`crates/spatial-stata/c/stplugin.c`, `stplugin.h`, `stata_bridge.c`) and persistent loaded `_scpc_rust_plugin` (`ping` + `exact_scpc`) with two-pass Modified Gram-Schmidt (`crates/spatial-stata/src/lib.rs`) on `[X, 1]` / `[X]` to build `proj_u`/`proj_v` inside Rust.
   - `scpc-stata/src/scpc.ado` adds `backend(string)` (honoring `$SPATIAL_BACKEND`) and `nthreads(integer -1)`, plus `ereturn local backend`. Requires `e(cmd) == "regress"`, unclustered, unweighted, numeric non-omitted/non-factor regressors, `5 <= N < 4500`, `0.001 < avc < 0.99`.

## 3. Changed and Added Files by Repository

- `spatial-rust` (`123dfcb`): `.gitignore`, `Cargo.toml`, `Cargo.lock`, `crates/spatial-core/{Cargo.toml,src/{basis,calibration,covariance,distance,lib,linalg,quadrature,scpc_op,student_t,threads}.rs}`, `crates/spatial-py/{Cargo.toml,pyproject.toml,src/lib.rs}`, `crates/spatial-r/{Cargo.toml,src/lib.rs}`, `crates/spatial-stata/{Cargo.toml,Makefile,c/{stata_bridge.c,stplugin.c,stplugin.h},src/{lib.rs,stata_spi.rs}}`, `fixtures/{generate_m1_fixtures.py,m1_exact_scpc.json,benchmark_m1.py,benchmark_m1_results.json}`.
- `scpc-python` (`9c9e062`): `pyproject.toml`, `src/scpc/core.py`, `src/scpc/types.py`, `tests/test_rust_m1.py`.
- `scpcR` (`cb35c73`): `.Rbuildignore`, `.gitignore`, `DESCRIPTION`, `NAMESPACE`, `R/core.R`, `man/scpc.Rd`, `src/{Makevars,Makevars.win,init.c}`, `src/rust/{Cargo.toml,Cargo.lock,src/lib.rs}`, `tests/testthat/test-rust-m1.R`.
- `scpc-stata` (`c1108e7` stacked on `9a84d46`): `src/scpc.ado`, `src/scpc.pkg`, `src/scpc_rust.plugin`, `tests/test_rust_m1.do`.

## 4. Verification Commands & Observed Results

All checks below were re-run cleanly after committing:

### Rust Workspace (`spatial-rust`)
```bash
cargo test --manifest-path spatial-rust/Cargo.toml --offline -p spatial-core --release
# 4 passed; 0 failed (0.07s)

cargo test --manifest-path spatial-rust/Cargo.toml --offline -p spatial-core --release --features portable
# 4 passed; 0 failed (0.24s)

cargo test --manifest-path spatial-rust/Cargo.toml --offline -p spatial-r
# 2 passed; 0 failed (0.00s)

cargo test --manifest-path spatial-rust/Cargo.toml --offline -p spatial-stata --locked
# 3 passed; 0 failed (0.09s)
```

### Python (`scpc-python`)
```bash
scpc-python/.venv/bin/pytest scpc-python/tests/test_rust_m1.py -q
# 21 passed in 7.63s (max |diff| vs pinned m1_exact_scpc.json and fallback <= 1e-10 asserted)

scpc-python/.venv/bin/ruff check scpc-python/src scpc-python/tests/test_rust_m1.py
# All checks passed!

scpc-python/.venv/bin/ty check --python scpc-python/.venv scpc-python/src
# All checks passed!

env PATH="/usr/bin:/bin:/usr/sbin:/sbin" scpc-python/.venv/bin/pytest scpc-python/tests -q
env PATH="/usr/bin:/bin:/usr/sbin:/sbin" SPATIAL_BACKEND=fallback scpc-python/.venv/bin/pytest scpc-python/tests -q
# Both runs: 108 passed, 57 skipped, 3 warnings (the 57 skipped tests require Rscript / fixest parity)
```

### R (`scpcR`)
```bash
R_LIBS_USER=.Rlib Rscript -e 'res <- testthat::test_dir("scpcR/tests/testthat", filter = "rust-m1", env = new.env(parent = asNamespace("scpcR")), reporter = "silent"); df <- as.data.frame(res); print(colSums(df[, c("nb", "failed", "skipped", "error", "warning")]))'
# nb=137, failed=0, skipped=0, error=0, warning=0

otool -L .Rlib/scpcR/libs/scpcR.so
# Only Accelerate, libSystem.B.dylib, libR.dylib, CoreFoundation (zero C++ stdlib linkage)
```
Full R `testthat` suite against `.Rlib` (`R_LIBS_USER=.Rlib Rscript -e 'testthat::test_dir("scpcR/tests/testthat", env = new.env(parent = asNamespace("scpcR")))'`):
- `rust-m1` and all non-`fixest` OLS unit tests pass.
- 6 skipped tests (1 `SCPCR_RUN_HEAVY_TESTS`, 5 `skip_on_cran` Stata parity).
- 17 pre-existing failures (`10+` cap reported) solely in `fixest` IV tests (`test_get_fixest_iv_design`, `test_get_scpc_model_matrix`, `test_is_fixest_iv_second_stage`, `test_scpc_fixest_iv_smoke`, `test_scpc`) because this machine has `fixest 0.12.1` installed while `DESCRIPTION` specifies `fixest (>= 0.14.0)`.

### Stata (`scpc-stata`)
```bash
cd scpc-stata/tests && /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b do test_rust_m1.do
# Exit code 0: "ALL MILESTONE 1 STATA TESTS PASSED!"
# Observed parity across bench_200 & bench_1000 (Euclid & latlong, cond & uncond, cvs):
#   stats max diff vs fallback <= 7.82e-14; cvs max diff <= 3.55e-15
#   bench_200 double precision vs pinned M1 JSON: 2.49e-14
#   nthreads(1) vs nthreads(4): stats diff = 0.00e+00, cvs diff = 0.00e+00
```

## 5. Benchmarks (`spatial-rust/fixtures/benchmark_m1_results.json`)

Workload: `y ~ x1 + x2`, Euclidean (`s_1`, `s_2`), `avc = 0.03`, `uncond = False`, `cvs = False`. Warm = median of 5 calls after 1 cold call.

| Dataset | Host | `rust_4t` Cold (ms) | `rust_4t` Warm (ms) | `rust_1t` Cold (ms) | `rust_1t` Warm (ms) | `fallback` Cold (ms) | `fallback` Warm (ms) | `q` |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `bench_200` | Python | 15.87 | 10.90 | 26.22 | 23.50 | 131.05 | 125.41 | 19 |
| `bench_200` | R | 38.00 | 11.00 | 42.00 | 24.00 | 463.00 | 131.00 | 19 |
| `bench_200` | Stata | 22.00 | 12.00 | 32.00 | 25.00 | 212.00 | 188.00 | 19 |
| `bench_1000` | Python | 56.32 | 51.97 | 160.51 | 137.96 | 576.35 | 572.24 | 8 |
| `bench_1000` | R | 73.00 | 53.00 | 163.00 | 138.00 | 1457.00 | 1146.00 | 8 |
| `bench_1000` | Stata | 67.00 | 53.00 | 144.00 | 136.00 | 1149.00 | 1120.00 | 8 |

For reference against prior `dev-cpp` warm medians recorded in `verification/reports/shared_rust_core_benchmarks.md` (`N=200 / 1000`: Python `38.9 / 277.6 ms`, R `16 / 61 ms`, Stata `105 / 1,227 ms`), Milestone 1 Rust beats all three baselines on both sizes.

## 6. Known Caveats & Next Scope (Milestone 2+)

- **Platform coverage**: Verified on macOS arm64 (`Darwin`). Linux and Windows native artifacts (`scpcR` on Linux/Windows, `scpc-stata` `.plugin` cross-builds, non-Apple wheels) are not yet tested natively here (only `--features portable` unit tests on macOS arm64).
- **`scpcR` standalone path dependency**: `scpcR/src/rust/Cargo.toml` uses relative workspace paths (`../../../spatial-rust/crates/...`), which works in the monorepo checkout and will need vendoring/packaging adjustments for standalone `R CMD build` / CRAN tarballs.
- **Host R parity environment**: Python R-parity tests (57 skipped without `Rscript`) and R's 17 `fixest` IV tests require an updated `fixest >= 0.14.0` + working `stringmagic` installation on the host.
- **Deferred beyond Milestone 1**: SCPC clustered / absorbed FE (`fixef_id`) / IV models / Nyström approximate large-N (`method = "approx"`), and all SPUR functionality (`spur-python`, `spuR`, `spur-stata` on `rust-core`).
