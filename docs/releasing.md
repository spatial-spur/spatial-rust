# Build and validate a matched set

Build from a single source revision. The version in the Python and R
distribution metadata identifies the frontends. The Rust core reports the
methods version and provider. Record release filenames, SHA-256 hashes,
and the revision with the manifest script.

```sh
python scripts/release_manifest.py dist/*.whl dist/*.tar.gz dist/r/*.tar.gz
```

The ordinary Python wheel installs the native extension with the frontend.
Build Linux wheels using a manylinux build image and audit the wheel tag and
external dependencies. Build the Windows and macOS wheels on their named
targets. The source archive includes the entire locked Cargo workspace.
Verify that it can build a wheel with the lock enforced.

The R archive uses a two-member workspace with the same dependency versions.
It filters and vendors locked registry sources with `cargo-vendor-filterer`,
compresses them, and tests the result with an empty Cargo cache and network
disabled. The archive must stay below the 5 MB CRAN source budget.

Install test artifacts from a working directory outside all source
checkouts. Use SCPC and SPUR with `backend="rust"` or `backend(rust)`.
The Python suite lives under `tests/python`. R installs run
`r/spatialSpur/tests/installed.R`. Licensed Stata runs use
`tests/stata/installed.do PACKAGE_DIRECTORY EMPTY_PLUS_DIRECTORY`.
The Stata executable must load the plugin and pass those tests before a
bundle is advertised for its target.

Run shared draw, Chetty, and historical verification after installing the
three artifacts. The Python process uses the installed package, and the R
process inherits its installed library.

```sh
python tests/parity/run.py --work "/tmp/spatial-spur-parity" \
  --stata-package "$PWD/dist/stata" --chetty
```

The optional `--historical-reference` flag compares coefficients and standard
errors with the archived authors’ MATLAB output and reports differences in
interval endpoints and decisions. The analytic infinite-decay endpoint can
change the reported cutoffs relative to finite-grid author arithmetic.

Verification outside the installed hosts also uses the Rust provider suites
and the independent small-matrix NumPy/SciPy reference. Compare exact and
Nyström results with recorded landmarks and seeds, especially on designs
with sparse regressors. Use representative local timing and peak memory
measurements for exact SCPC, approximate SCPC, local transforms, the
combined workflow, and half-life.
