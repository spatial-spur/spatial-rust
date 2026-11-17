# SPUR half-life of spatial persistence

Methods specification 0.1.0. This is the scale-invariant confidence set for the half-life of a scalar spatial local-to-unity process. The mathematical target uses the Lévy–Brownian basis of Müller and Watson and the population inversion in their equation (20). Source-specific approximation and reporting rules are recorded in [Conventions and decisions](conventions.md#spur-half-life).

## Geometry, basis, and null

For observed locations $s_i$, use Euclidean distances $D_{ij}=\lVert s_i-s_j\rVert$, diameter $d_{\max}=\max_{ij}D_{ij}>0$, and $M_1=I-\mathbf1\mathbf1^\top/n$. The Lévy–Brownian covariance at these locations has entries $(\lVert s_i\rVert+\lVert s_j\rVert-D_{ij})/2$, so its centered covariance is

$$
B_L=-\tfrac12 M_1DM_1.
$$

Let $R_q\in\mathbb R^{n\times q}$ have unit-norm orthogonal columns spanning the leading $q$ centered Lévy–Brownian eigendirections. The paper scales these columns by $\sqrt n$. The choice of unit columns rescales the observation and projected covariances consistently. Set $z=R_q^\top y=R_q^\top M_1y$.

For $h>0$, $c=\log(2)/h$. Under the canonical Gaussian null at $h$, the projected covariance, up to a common positive variance scale, is

$$
\Omega(h)=R_q^\top\frac{\exp(-cD)}{2c}R_q,
$$

where the exponential is entrywise. The regular experiment uses $q$ with positive retained covariance directions, positive-definite $\Omega(h)$ at each tested or mixed $h$, and $z\ne0$. The limiting projected covariance as $h\to\infty$ is $-R_q^\top D R_q/2$, by expansion of the exponential and $R_q^\top\mathbf1=0$. As $h\to 0$, the covariance *shape* approaches the identity at distinct sites. With duplicates the limiting shape instead comes from the projected co-location indicator. Boundary tests require a well-defined limiting law. [Müller and Watson 2024, §2.2 equation (2), §§2.5 and 4.1, and §4.4 equation (20)](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf)

## Angular likelihood and uniform alternative

For nonzero $z$ and positive-definite $\Omega$, write $u=z/\lVert z\rVert$. Up to a constant depending only on $q$, the Gaussian direction has density

$$
f_\Omega(u)=\lvert\Omega\rvert^{-1/2}(u^\top\Omega^{-1}u)^{-q/2}.
$$

Multiplying $\Omega$ by a positive scalar leaves this density unchanged. In particular, dropping $1/(2c)$ at each positive $c$ leaves the angular comparison unchanged when both parts of the formula use the same covariance shape. The alternative weights half-lives uniformly from zero to the sample diameter. A normalized score for testing a particular $h_0$ is

$$
L_{h_0}(u)=
\frac{\frac{1}{d_{\max}}\int_0^{d_{\max}} f_{\Omega(h)}(u)\,dh}
     {f_{\Omega(h_0)}(u)}.
$$

Equation (20) uses the integral without its uniform normalizer and omits the null determinant factor. The two scores differ by a positive constant for each $h_0$, hence have the same ranking and yield the same test when calibrated under that null. The integral averages densities, and it uses a uniform measure on $h$. The tested null $h_0$ may lie beyond $d_{\max}$. [Müller and Watson 2024, §4.4, equation (20)](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf#page=16)

## Confidence-set inversion

For each null $h_0$, the theoretical null direction is distributed as $Z_0/\lVert Z_0\rVert$ for $Z_0\sim N(0,\Omega(h_0))$. The critical value is the population $(1-\alpha)$-quantile of $L_{h_0}$ under that null. Reject when $L_{h_0}(u)$ exceeds the critical value. The confidence set contains all null half-lives whose tests do not reject. When the null score CDF is continuous and strictly increasing at its critical value, this includes the corresponding p-value equal to $\alpha$.

The confidence set can be disconnected. Taking the hull from its lowest to highest accepted value is a reporting rule. Treat a tested or reported infinite endpoint, an empty accepted set, and an undefined likelihood as distinct states. The population quantile rule is defined even when all likelihood ratios tie. With $q=1$, or when all retained covariance shapes are proportional, the direction carries no information about the covariance shape. A finite simulation can count no strict exceedances in that case, so the finite p-value rule needs a separately documented convention. [Müller and Watson 2024, §§4.1 and 4.4](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_SpatialUnitRoots.pdf#page=14)

## Independent checks

Adding a constant to $y$ leaves $z$ unchanged. Multiplying $y$ by a nonzero scalar leaves the score unchanged, since it is even in direction. An orthogonal change of coordinates within the same retained subspace preserves the score when applied to the observation and all null and alternative covariances. Changing Euclidean distance units by a positive factor $a$ multiplies half-lives and the diameter by $a$, with unchanged normalized inference.

A three-site, two-direction example has a closed-form confidence endpoint and an independently integrated angular likelihood. Its normalized population upper endpoint exceeds 3, making it useful for distinguishing the continuous method from a grid with a gap between 3 and an open-endpoint sentinel. [Reference calculations](reference/README.md)
