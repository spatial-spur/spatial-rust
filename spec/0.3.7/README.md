# SCPC and SPUR methods specification 0.3.7

This record extends [0.3.6](../0.3.6/README.md). The statistical target for half-life is the [continuous confidence set](../0.1.0/spur-half-life.md) in Müller and Watson (2024, equation 20).

## SPUR covariance and components

For a levels or residual projection $M$ that removes the constant direction, the Lévy–Brownian covariance is $-MDM/2$, where $D$ contains distances divided by sample diameter. The residual basis, the residual I(1) null, and the half-life basis use this covariance directly. Residual-response projection follows 0.3.6.

For fixed SPUR dimension $q$, require a positive retained covariance rank. If the gap between eigenvalues $q$ and $q+1$ is at most $1.6\times10^{-11}\max(1,\lambda_1)$ in the normalized covariance, report an ambiguous component count. A basis may contain a full tied eigenspace. Its mathematical likelihood is invariant to rotations within that subspace. Simulating with a fixed seed in different bases can produce different finite p-values.

## SPUR simulation and critical values

Use $B$ independently generated standard normal columns for the null comparison. The I(0) and I(1) statistics also use an independent block of $B$ columns to choose the alternative. The native request stores null draws in the first $q$ rows and tuning draws in the next $q$ rows. Methods that share a sample can reuse each block across tests. If the tuning search cannot bracket 50% estimated power on the high side, use the strongest trial according to the tuning simulation.

For an upper-tail comparison at a specified null, compute $p=(1+\#\{T_b\ge T_{\rm obs}\})/(B+1)$. The $\alpha$ critical value is the $\lfloor\alpha(B+1)\rfloor$-th largest simulated statistic. Reject above that cutoff. A zero rank represents an unreachable rejection level and is returned as infinity, or missing in Stata. A valid rank test assumes exchangeability under the specified null, conditional on the geometry, eligible design, and independent tuning draws.

The I(0) calculation takes maxima over its finite nuisance candidates, including the exact weak-correlation limit at infinite decay with covariance kernel $1[D_{ij}=0]$. The interval between nuisance candidates still contributes to the theoretical supremum.

## Half-life inversion

Use half-life $h$ in sample-diameter units and search on $t=h/(1+h)$. Evaluate the two boundary laws explicitly. At $h=0$, the covariance kernel is $1[D_{ij}=0]$. At $h=\infty$, use the projected Lévy–Brownian covariance. Recolor the same null draws at each candidate. The score averages the angular Gaussian alternative densities using five-point Gauss–Legendre integration on each of the intervals with edges $0,2^{-9},2^{-8},\ldots,1$.

Probe $t=i/128$, $i=0,\ldots,128$, and bisect each detected acceptance transition to a bracket width of $10^{-7}$ in $t$. Round detected components outward within their brackets. The search has finite resolution. An acceptance component or rejection gap between all probes can escape detection, and the $t$-scale tolerance maps to wider half-life brackets near infinity.

Return the detected components and their hull. A set can be empty, bounded, or unbounded. Report infinity for an unbounded endpoint, with a separate indicator in Stata. The initial probes and their p-values are available in the Python and R results.

The shared adapters use SPUR API 5 for the two draw blocks and the half-life result. SCPC uses API 7.

## Calibration grid audit

On a six-site conditional OLS fixture, SCPC selects three components and returns a 95% cutoff of 29.762524. Independent Gaussian integration gives a rejection probability of 5.006689% at decay 0.803351. That decay lies between the first two points on the geometric grid. The fixture demonstrates an interior omission even with the finite upper endpoint included. The SCPC calibration still uses a finite grid and stops at a small positive correlation above the co-location floor.

For a fourteen-site SPUR residual fixture with two retained directions, the legacy I(0) grid gives 5.012716% rejection at independence when calibrated to 5% on the grid using population probabilities. Including the exact weak-correlation limit reduces rejection in that case to 5%. The test still samples finitely many interior nuisance values. The positive grid starts at 0.0001, so the current implementation requires a co-location floor below that value.

Run [the grid audit](audit_inference_grids.py) with NumPy, SciPy, and the matching development checkouts of `scpc-python` and `spur-python` on the Python path.
