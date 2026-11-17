# Shared-core implementation conventions

Methods specification 0.2.0. These choices describe the corrected package profile, separately from the inherited population theory. They supersede the corresponding historical package descriptions in version 0.1.0.

## SCPC

- Every native and fallback OLS path uses the version 0.1.0 required projection order. IV uses the transposed residual maker and weight order in [Conditional IV SCPC](scpc-iv.md). Absorbed effects belong to the full observation-level residual projection before aggregation.
- Critical values use the upper bisection bracket endpoint. Conditional intervals use the maximum of the location and conditional cutoffs, and p-values use the maximum of the two family p-values. The existing finite decay grid and component-selection profile remain in use.
- Automatic dispatch may use a host fallback for an unsupported model. An explicitly requested Rust backend must report unsupported inputs. The backend choice and thread count are separate options. One thread remains a Rust computation.

## SPUR validity and finite simulations

- Require at least two retained directions for SPUR inference. With one direction the angular likelihood has no covariance-shape information. Report an input/degeneracy error rather than a finite interval or simulation p-value. Require residual tests to have `q <= n-rank(X)` and require the retained basis covariance to have positive numerical rank.
- Reject wholly proportional projected covariance shapes as an uninformative experiment. The half-life implementation compares trace-normalized shapes across its actual null and alternative grids. I(0) compares its two mixture components. The I(1) diagnostic compares shapes at average correlations 0.01, 0.5, and 0.95. This three-point diagnostic is a numerical safeguard, not a proof of identifiability for every geometry. The normalized-shape tolerance is `1e-12`.
- Normalize nonzero design columns before rank checks, use two orthogonalization passes, and apply a dimensionless `1e-10` residual-norm threshold. Check positive retained basis eigenvalues relative to the leading value (`1e-12`). Reject constant responses and zero projected responses. Center and scale responses before likelihood calculations, so changing a response's level or nonzero units does not create a degeneracy error.
- Keep the author's 131 null / 50 alternative half-life grid, common Gaussian draw matrix, strict simulation exceedance p-value, and strict grid acceptance rule. An empty accepted set, an accepted open-upper sentinel, and an invalid experiment have distinct internal outcomes. Public missing-value representations retain their package conventions.
- Canonicalize each isolated eigenvector by making its earliest largest-magnitude entry positive. The host generates the same Gaussian matrix before either backend, and fallback eigensolvers must not consume additional host RNG draws. For a fixed host seed and isolated eigenvalues, compare public native/fallback results and final host RNG state. Seeds are host-specific, so an equal integer seed across languages is not a common-draw test.
- Eigenvector signs do not identify a basis inside a repeated eigenspace. Native operations report the tied-eigenspace flag and automatic dispatch preserves it when selecting a host fallback. Forced host fallback calls do not independently diagnose ties. Exact finite-simulation parity across independently computed tied bases is not promised. Compare invariant calculations or supply a common basis when studying such cases.
- Limit each diagnostic power-calibration bracketing phase to 200 evaluations, with a calibration error if a bracket cannot be found. A legal experiment can fail to reach the selected 50% power target, which is distinct from invalid rank or proportional covariance shapes.
- Use type-7 linear empirical quantiles in all three implementations. Stata's earlier discrete-index fallback quantile has been replaced to match its native calculation. This and the corrected basis orientation can change a finite seeded realization without changing the target Gaussian draw law.

## Distance units

A normalized half-life is a fraction of the sample diameter. For Euclidean coordinates, the reported diameter and absolute endpoints use those coordinates' units. For geographic inputs, the existing angular distance convention is converted to metres exactly once using the package Earth-circumference constant. Scaling Euclidean coordinates by a positive factor scales absolute half-lives and the diameter by that factor, while leaving normalized inference unchanged. The historical unconditional Earth conversion in the R and Stata Euclidean wrappers is a corrected defect.

The error policies above are adopted for the clean shared-core branches. They do not rewrite archived release behavior or change the population method definitions.

The corrected independent Python extension requires `SPUR_API_VERSION = 2`; Stata requires `ping_spur_v2`. Automatic dispatch routes older native artifacts to the corrected host fallback. An explicit Rust request with an older artifact reports an incompatibility.
