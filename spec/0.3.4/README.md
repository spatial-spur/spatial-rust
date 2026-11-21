# SCPC and SPUR methods specification 0.3.4

This record extends the package conventions in [0.3.3](../0.3.3/README.md). The conditional SCPC interval, Gaussian rejection rule, and SPUR diagnostic and half-life targets use the definitions in the earlier versions.

## SPUR transformations

For nearest-neighbor differencing, average the observations at the smallest positive distance, with equal weights across ties. A sample with only one distinct location cannot use this transformation. The isotropic transformation uses all neighbors with $0 < d \leq b$ for radius $b$. A row with no such neighbors has a transformed value of zero. These choices adopt the R and Stata conventions for Python.

Measure geographic radii and half-life endpoints in meters with Earth radius $6{,}371{,}000.009$ m. Use $2\pi R$ to convert fractions of the circumference. Normalized distances continue to use the sample diameter.

Apply the LBM GLS operator to the evaluated model columns and remove any residual column means from the transformed values. Both Python and R use the fitted design, including its intercept choice, for residual diagnostics and refit transformed no-intercept formulas without a constant. Stata's SPUR regression command fits with an intercept. Separate SPUR diagnostic and transformation commands can be used on explicit variables in Stata.

## Development bindings

The native host bindings use SCPC model API 7. It takes fitted coefficients, observation scores, the corresponding bread, and the selected design and site information needed for inference. Python exports `run_model_scpc`, R registers `wrap__rust_run_model_scpc`, and Stata dispatches `model_exact_scpc` and `model_approx_scpc`. The Rust crate also exposes this operation and its common `ScpcOutput` at the crate root. The host packages provide pure-language inference and their usual backend selection.
