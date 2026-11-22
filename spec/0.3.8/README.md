# Methods specification 0.3.8

This record updates the numerical conventions in [0.3.7](../0.3.7/README.md). The statistical targets remain those in [the theory specification](../README.md).

## Half-life likelihood and inversion

The mixture integrates the angular Gaussian density uniformly over half-lives from zero to one sample diameter. Let $d_{\min}$ be the smallest strictly positive normalized inter-site distance. The first Gauss–Legendre panel starts at zero and ends at the largest power of two at or below $\min(2^{-9},d_{\min}/32)$, subject to the normal binary64 floor. Subsequent panels double that endpoint up to one. Each panel uses five quadrature nodes. The partition depends only on geometry and applies to the observed and simulated scores.

The inversion evaluates both endpoint laws and samples $t=h/(1+h)$ in steps of $1/128$. Below the first positive probe it adds four probes per halving of $h$ until the probes reach $d_{\min}/32$, subject to the same floor. Brackets around changes in acceptance are refined to the smaller of $10^{-7}$ in $t$ and $10^{-5}$ of their initial width. Detected components are rounded outward at their respective brackets. The package reports the components and their hull as a numerical approximation. Narrow components or gaps between probes can escape detection. Python and R also return the initial probe locations and rank p-values.

For sites $(0,10^{-8},0.3,0.6,1)$ and the response contrast $(1,-1,0,0,0)$ in the full centered space, independent adaptive quadrature gives mixture density $3.63291497\times10^{-9}$ with the common angular constant omitted. The revised quadrature matches it within relative tolerance $10^{-5}$. An eight-site regression fixture exercises a short accepted component of order $10^{-6}$ between zero and the first positive probe in version 0.3.7.

## SPUR diagnostics and transformations

- For residual tests, $M_X=I-P_X$ and the design must span the constant to the existing numerical tolerance. For an offset $a$, prepare $y-aP_X\mathbf{1}=(y-a\mathbf{1})+aM_X\mathbf{1}$. A residual of the unit constant below eight machine epsilons is treated as roundoff. Otherwise the correction is carried into projection and the response variation check allows for its roundoff. The levels tests and cluster demeaning also subtract a safe offset before computing means. The intercept-plus-indicators and full-indicators representations agree on the large-offset regression fixture.
- Retaining $q$ directions requires a gap larger than $\max(1.6\times10^{-11}\lambda_1,8n\epsilon)$ at the cutoff in the normalized LBM covariance. The term $8n\epsilon$ allows for construction roundoff on the scale of the unprojected kernel. The relative term follows the projected spectral scale. A cutoff inside a numerically tied eigenspace remains ambiguous.
- If the co-located-pair fraction is $f$, the I(0) benchmark uses the attainable average correlations from $f$ through $0.03$. The endpoint covariance at $f$ is $\mathbf{1}[D=0]$. The interior grid starts above $f$ and includes $0.03$ whenever $f<0.03$. I(0) uses the limiting covariance for tuning when the $0.001$ reference is unattainable. I(1) places its correlation anchors at fractions of the attainable interval $(f,1)$. When a tuning search reaches a kernel whose positive-distance correlations have rounded to zero, doubling the decay further leaves that covariance unchanged. The search can return the strongest tested alternative if it has yet to reach its target.
- The I(0) composite calibration takes maxima on its numerical grid. The point-null Monte Carlo rank guarantee and the finite-grid approximation to composite calibration have different scopes.
- Cluster demeaning uses only group membership and accepts a one-observation sample. SPUR and half-life diagnostics require at least two retained directions and sufficient residual dimension. LBM-GLS retains normalized covariance eigenvalues greater than $10^{-10}$.

## SCPC and interpretation

- Python reconstructs PyFixest columns with the fitted Formulaic transformation state and retained sample. Fixed-effect OLS and IV tests compare the result to a fit on the corresponding precomputed transformed columns.
- Euclidean Nyström SCPC subtracts a finite location offset and scales before forming the centered coordinate cross-product. Its landmark consolidation still selects a leading direction at each step.
- Exact SCPC computes without sampling the spatial matrix. Both SCPC methods approximate the supremum over the decay parameter on a finite grid. The [independent grid audit](../0.3.7/audit_inference_grids.py) records an off-grid rejection probability of 5.006689% at a nominal 5% on one six-site geometry.
- The prespecified levels and transformed analyses each have their own inference assumptions. Selecting a branch using a SPUR diagnostic on the same observations requires a separate coverage argument.

The native adapters continue to use SCPC API 7 and SPUR API 5. The host packages require methods version 0.3.8 for native parity with these conventions.
