# Exact conditional OLS SCPC

Methods specification 0.1.0. One coefficient, one observation per location index, and a full-rank, unweighted OLS design. Distances in the published setup are Euclidean. Exact refers to the spatial covariance and eigenspace on all sample locations. Conditional inference here concerns the Gaussian benchmark conditional on the observed design and locations. See the [source register](sources.md) for the unconditional SCPC source and the regression extension.

## Coefficient and scores

Let $y\in\mathbb R^n$ be the outcome and $V=[x_{\mathrm{raw}},Z]$ the design, with $Z$ containing the other regressors, including an intercept if fitted. Write $M_A$ for orthogonal projection onto the complement of the column space of $A$, and $M_A=I$ if $A$ has no columns. Use the Frisch–Waugh–Lovell regressor $x=M_Zx_{\mathrm{raw}}$, so $V$ and $[x,Z]$ span the same space. For a hypothesized value $\beta_0$,

$$
\widehat\beta=\frac{x^\top y}{x^\top x},\qquad
\widehat e=M_Vy,\qquad
\widehat u=\operatorname{diag}(x)\widehat e,\qquad
S_{xx}=\frac{x^\top x}{n}.
$$

The benchmark true score is $u=\operatorname{diag}(x)e$. In particular $n^{-1/2}\mathbf 1^\top u=S_{xx}\sqrt n(\widehat\beta-\beta_0)$ under the null in the OLS model. The variance estimator below uses the fitted residuals. [Müller and Watson 2023, §§2.1–2.2, equations (3) and (5)](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=4)

## Spatial basis and statistic

For $D_{ij}=\lVert s_i-s_j\rVert$, set $K(c)_{ij}=\exp(-cD_{ij})$ and

$$
\bar\rho(c)=\frac{1}{n(n-1)}\sum_{i\ne j}K(c)_{ij},\qquad
\bar\rho(c_{\min})=\bar\rho_{\max},\qquad
M_1=I-\frac{\mathbf 1\mathbf 1^\top}{n}.
$$

Here $\bar\rho_{\max}$ is the user's correlation bound. With at least one positive off-diagonal distance, $\bar\rho(c)$ decreases strictly from 1 to the off-diagonal fraction of zero distances. A finite positive $c_{\min}$ exists for a bound strictly between those limits. At $c\to\infty$, $K(c)\to I$ for distinct sites and tends to the indicator of co-location for repeated sites. The i.i.d. model used to select $q$ is separately defined by $I$.

Let $B_q$ contain the leading $q$ unit-norm orthogonal eigenvectors of $M_1K(c_{\min})M_1$ in the centered subspace. The papers write the corresponding eigenvectors as $R=\sqrt n B_q$. In unit-norm notation the same quantities are

$$
\widehat\sigma_q^2=\frac{\lVert B_q^\top\widehat u\rVert^2}{q},\qquad
SE_q=\frac{\widehat\sigma_q}{S_{xx}\sqrt n},\qquad
t_q(\beta_0)=\frac{\widehat\beta-\beta_0}{SE_q}.
$$

These formulas use $q$ nonconstant directions. Centering $\widehat u$ leaves the projections unchanged. The regular case has a full-rank design, an identifiable coefficient, and a positive projected variance for the observed statistic. If a selected prefix cuts through a tied eigenvalue, the prefix subspace may be ambiguous. A joint orthogonal change of basis within a fully retained subspace preserves the statistic and the Gaussian rejection probability at fixed $q$. [Müller and Watson 2022, §2, equations (3)–(6) and Remark 2.3](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf#page=5) · [Müller and Watson 2023, §2.2](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=5)

## Two Gaussian benchmark families

Set $b_0=\mathbf 1/\sqrt n$, $D_x=\operatorname{diag}(x)$, and $S_x=\operatorname{diag}(\operatorname{sign}(x))$. The following expressions apply directly when the entries of $x$ are nonzero. The convention at an exact zero is explicitly tracked in the [convention record](conventions.md#scpc).

| Family | Model for the numerator and denominator | Weight matrix |
| --- | --- | --- |
| Location benchmark | A Gaussian score vector with covariance $K(c)$. Centering for the denominator is absorbed by $B_q^\top\mathbf 1=0$. | $W_{\mathrm{loc}}=[b_0,\ B_q]$ |
| Conditional OLS benchmark | $a\mid(V,s)\sim N(0,K(c))$, $e=S_xa$, $u=D_xe$, and $\widehat u=D_xM_Ve$. | $W_{\mathrm{cond}}=[S_xD_xb_0,\ S_xM_VD_xB_q]$ |

In the conditional model the numerator is $b_0^\top D_x S_xa$. The denominator projections are $B_q^\top D_xM_VS_xa$. Transposing these row maps gives the columns shown above. [Müller and Watson 2023, §2.3, equation (7), and §5.1 following equation (8)](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=14)

For either family $j$, let $G\sim N(0,\Omega_j(c))$ with $\Omega_j(c)=W_j^\top K(c)W_j$. At a cutoff $k\ge 0$ in ordinary t units, the null rejection event is

$$
G_0^2-\frac{k^2}{q}\lVert G_{1:q}\rVert^2>0.
$$

This is equation (8) on substituting the unit-norm weights. Equivalently, it is the sign of the Gaussian quadratic form with symmetric matrix $\Omega_j(c)^{1/2}\operatorname{diag}(1,-k^2 I_q/q)\Omega_j(c)^{1/2}$. The event itself also defines the probability for a singular Gaussian covariance. It avoids dividing by a simulated denominator that might be zero.

Define $R_j(k)=\sup_{c\ge c_{\min}}\Pr_{j,c}(\text{reject at }k)$ and $R(k)=\max_j R_j(k)$. The conditional SCPC critical value is the smallest size-controlling cutoff, expressed as

$$
k_*=\inf\{k\ge 0:R(k)\le\alpha\}.
$$

In continuous regular cases it equals the larger of the two family critical values. Report $[\widehat\beta-k_*SE_q,\ \widehat\beta+k_*SE_q]$. A two-sided p-value for the same combined test is $R(|t_q(\beta_0)|)$ in the continuous case. Retain family p-values separately if reporting them. If the Gaussian law is degenerate, use the stated strict rejection event and make any p-value tie convention explicit.

This specification uses the continuum of null decays. A finite set of decays yields size control over that set. The mathematical critical value does not prescribe a bisection step, a quadrature rule, or a finite grid. [Müller and Watson 2023, §§2.3 and 5.1](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=6) · [Müller and Watson 2022, equation (5)](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf#page=6)

## Required projection order

Every package and backend implementing conditional OLS SCPC must use the denominator weight columns $S_xM_VD_xB_q$. Multiply the basis columns by the residualized regressor, apply the residual maker for the full fitted OLS design, then multiply by the signs of that regressor. If the original fit uses no intercept, $M_V$ uses that same design. This requirement applies to the shared Rust core and to Python, R, and Stata implementations, including fallbacks. The [fixed-component reference calculation](reference/README.md#conditional-scpc-with-a-fixed-component-count) distinguishes the two orders.

For a clustered OLS extension, let $\Xi$ contain one column per cluster with the residualized regressor values for that cluster and zeros elsewhere. For clusters with nonzero regressor norm, $\Xi_s$ normalizes each column to unit length. Use $\Xi_s^\top M_V\Xi B_q$ with a cluster-level basis $B_q$ and an observation-level residual maker $M_V$. This reduces to the unclustered order for one observation per cluster with nonzero $x$. The numerical specification and fixtures in version 0.1.0 cover unclustered OLS. [Müller and Watson 2023, §5.1, p. 1062](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=14) · [Clustered MATLAB source](sources.md#author-software-and-replication)

## Recommended choice of component count

For an automatically selected $q$, use the authors' recommended expected-length criterion for the 95% i.i.d. **location** benchmark. For each eligible $q$, let $k_{\mathrm{loc},q,.05}$ be the 5% location-family cutoff calibrated over $c\ge c_{\min}$. Minimize

$$
k_{\mathrm{loc},q,.05}\frac{\Gamma((q+1)/2)}{\sqrt q\,\Gamma(q/2)}
$$

over the admissible centered directions. The omitted multiplier is common to all $q$. The conditional critical value is evaluated after choosing $q$, and the same chosen $q$ serves other reported levels for the same geometry and correlation bound. A request with explicitly supplied $q$ uses that value. Document any finite candidate search, tie policy, or different selection level in its profile. [Müller and Watson 2022, equation (6)](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf#page=6) · [Müller and Watson 2023, §2.2](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=5)

## Independent checks

The identity covariance gives Student-$t_q$ for the location benchmark. Changing Euclidean distance units by a positive factor $a$ changes $c_{\min}$ by $1/a$, preserving the benchmark covariance. A four-site calculation checks the OLS statistic, calibrated decay, both population critical values, the combined interval, and the order of conditional weighting and projection. Its supplied $q=1$ isolates those operations from the automatic component search. [Reference calculations](reference/README.md)
