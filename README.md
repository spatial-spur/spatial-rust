# Spatial SPUR

Spatial SPUR provides spatial correlation-robust inference (SCPC), spatial unit root
diagnostics, half-life confidence sets, and spatial transformations. One Rust engine
implements methods version 0.4.0. Python, R, and Stata fit models in their own runtimes.

## Install from a matching build

| Frontend | Installation | Public entry points |
| --- | --- | --- |
| Python | install the `spatial_spur` wheel with `python -m pip install PATH_TO_WHEEL` | `from spatial_spur import scpc, spur` |
| R | `R CMD INSTALL spatialSpur_0.6.0.tar.gz` | `library(spatialSpur)` |
| Stata | `net install spatial_spur, from("/path/to/your/platform/bundle")` | `scpc`, `spur`, `spurtest`, `spurtransform`, `spurhalflife` |

Python wheels include the native library. Install the source checkout as
`pip install ".[pyfixest]"` if you use PyFixest. R installs need Cargo, Rust
1.85 or newer, and R’s normal compilation tools. The R release archive includes
locked Rust dependencies for offline compilation. Stata bundles contain one
platform-specific plugin shared by every command.

Build all three distributions from the repository root.

```sh
python -m pip install maturin
maturin build --release --locked
maturin sdist
scripts/build_r_source.sh dist/r
scripts/build_stata_bundle.sh dist/stata
```

The Stata link script supports macOS and Linux. The Windows linker is
`crates/spatial-stata/build_windows.ps1`. Use the Stata bundle tested for
your runtime and architecture. CI builds Python ABI3 wheels for CPython 3.11
and later on macOS, manylinux x86-64, and Windows x86-64. The runtime validation
in this branch uses Apple Silicon. See [release validation](docs/releasing.md).

## Fit once and supply the locations

### Python

```python
import statsmodels.formula.api as smf
from spatial_spur import load_chetty_data, scpc, spur

data = load_chetty_data().query("state not in ['AK', 'HI']")
data = data.dropna(subset=["am", "fracblack", "lat", "lon"])
fit = smf.ols("am ~ fracblack", data).fit()
print(scpc(fit, data, lon="lon", lat="lat", backend="rust"))
result = spur("am ~ fracblack", data, lon="lon", lat="lat", q=6, nrep=400)
print(result.summary())
```

### R

```r
library(spatialSpur)
d <- subset(load_chetty_data(), !(state %in% c("AK", "HI")))
d <- d[complete.cases(d[c("am", "fracblack", "lat", "lon")]), ]
fit <- lm(am ~ fracblack, d)
scpc(fit, d, lon = "lon", lat = "lat", backend = "rust")
p <- spur(am ~ fracblack, d, lon = "lon", lat = "lat", q = 6, nrep = 400)
p$transformed
```

### Stata

Given a dataset containing coordinates and regression variables, put the
latitude in `s_1` and the longitude in `s_2`.

```stata
generate double s_1 = lat
generate double s_2 = lon
regress am fracblack
scpc, latlong backend(rust) weights
spur am fracblack, latlong q(6) nrep(400)
```

The adapters use the fitted sample. For matching Python and R data frames,
preserve unique row identifiers. A PyFixest fit needs stored data and either
matching fitted row values or an explicit `row_id` column. For R fixest, use
`demeaned = TRUE, data.save = TRUE`. Explicit coordinate matrices are accepted
in fitted order. A spatial cluster must identify one actual site with constant
coordinates across its fitted observations.

SCPC returns coefficient-specific 95% critical values and separate
location-only `cv`. Python and R return spatial weights in the reported site
order. Stata’s `weights` option returns `e(W)` with `e(site_rows)`.
The full geometry and Nyström approximation are explicit choices. Approximate
results include the seed and landmark settings. An exact run has no
approximation seed.

The [current methods and support matrix](spec/0.4.0/README.md) describe the
statistical scope, model requirements, geometry, and randomization. The
[release instructions](docs/releasing.md) cover clean installations and
cross-language verification.
