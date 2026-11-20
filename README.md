# Shared native SCPC and SPUR build

This branch uses the stock OxiBLAS 0.2.4 crates for matrix multiplication and dense symmetric eigendecomposition on every target. The native core also uses Rust scalar exponentials and its existing PSD-aware weighted eigensolve for SCPC calibration. Rayon in the shared core controls the worker budget. The build uses OxiBLAS's default CPU dispatch, with its experimental SME2 feature off.

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

## Native symbols

The stock OxiBLAS crates supply CBLAS symbols. The Stata plugins export the two entry points in the Stata plugin ABI. The Python and R shared libraries retain the stock exports and use their host's normal extension loader. On Apple Silicon, loading these libraries preserved the host's existing BLAS symbol binding. This verification applies to the tested macOS artifacts and loaders.
