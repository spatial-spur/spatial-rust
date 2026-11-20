# Shared native SCPC and SPUR build

This branch uses the stock OxiBLAS 0.2.4 crates for matrix multiplication and dense symmetric eigendecomposition on every target. The native core also uses Rust scalar exponentials and its existing PSD-aware weighted eigensolve for SCPC calibration. Rayon in the shared core controls the worker budget. The default build uses OxiBLAS's default CPU dispatch.

The three adapters consume the same `spatial-core` sources and locked OxiBLAS revision. Build with Rust 1.89 or newer. The OxiBLAS license is reproduced from the [0.2.4 publication source](https://github.com/cool-japan/oxiblas/blob/995879190ba29ebf8264a4216e4d70e970a6da4f/LICENSE) in [LICENSES](LICENSES/OxiBLAS-Apache-2.0.txt) and included with each package's native distribution.

## Build and test

From this workspace:

```sh
cargo test --locked --release -p spatial-core -p spatial-stata
maturin build --locked --release --manifest-path crates/spatial-py/Cargo.toml
```

Install the resulting `spatial-spur-rust` wheel alongside either Python package. Both packages discover the shared extension through the existing Rust adapter.

Build each R source tarball from the corresponding `oxiblas-only` package checkout:

```sh
bash scripts/build_r_source.sh /path/to/scpcR /path/to/artifacts
bash scripts/build_r_source.sh /path/to/spuR /path/to/artifacts
```

The builder synchronizes the Rust sources and vendors the locked crates for offline installation. Set `SPATIAL_RUST=0` when installing a tarball to use the R implementation.

Build the Stata plugins from matching package checkouts. Set the Makefile's output paths explicitly when the checkouts are in separate worktree directories:

```sh
make -C crates/spatial-stata \
  PLUGIN_OUT=/path/to/scpc-stata/src/scpc_rust.plugin \
  SPUR_PLUGIN_OUT=/path/to/spur-stata/spur_rust.plugin
```

## Optional Arm SME2 build

SME2 is off by default. The opt-in builds use OxiBLAS's runtime checks to select SME2 kernels for eligible operations on supported Arm systems. Other hosts and operations use OxiBLAS's regular CPU dispatch.

Set `OXIBLAS_DISABLE_SME2=1` before starting the host process to use regular CPU dispatch with an opt-in build.

- **Python.** Add `--features sme2` to the Maturin command and install the resulting shared wheel.
- **R.** Build the source archives as above. Set `SPATIAL_SME2=1` when installing either archive with `R CMD INSTALL`. The same source archives support the default build.
- **Stata.** Add `SPATIAL_SME2=1` to the Make command above to build both plugins.

## Native symbols

The stock OxiBLAS crates supply CBLAS symbols. The Stata plugins export the two entry points in the Stata plugin ABI. The Python and R shared libraries retain the stock exports and use their host's normal extension loader. On Apple Silicon, loading these libraries preserved the host's existing BLAS symbol binding. This verification applies to the tested macOS artifacts and loaders.
