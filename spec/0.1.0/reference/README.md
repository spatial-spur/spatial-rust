# Independent reference calculations

Specification 0.1.0. These are small population-method cases with analytic derivations. They exercise exact arithmetic targets for inference without any release seed, library eigenspace convention at a split cutoff, or finite covariance grid. The historical conditional covariance in the first fixture is an intentional comparison against a separate formula.

Run the checker with Python 3.11 or newer, NumPy and SciPy. It reads the adjacent [fixtures](fixtures.json). It uses dense OLS, symmetric eigendecomposition, direct matrix products and independent numerical integration. It imports no SCPC, SPUR or Rust code. The tolerance in the fixture is for these checks in double precision.

```sh
# From the spatial-rust repository root
python spec/0.1.0/reference/check.py
```

## Conditional SCPC with a fixed component count

Take locations $0,1,2,3$, fit an intercept and $x=(-2,-1,1,2)^\top$, and use $y=(-1,-2,0,3)^\top$. Set $\bar\rho_{\max}=17/48$, $q=1$, $\beta_0=0$, and $\alpha=0.05$. Writing $r=\exp(-c)$, the average correlation is $(3r+2r^2+r^3)/6$. Thus $c_{\min}=\log 2$, with $0\le r\le 1/2$ in the benchmark including its limit. At $c_{\min}$, the centered eigenvalues are $1,\ 9/16,\ 3/8$, and the leading direction is $b=(2,1,-1,-2)^\top/\sqrt{10}$ up to sign. OLS gives $\widehat\beta=1$ and $\widehat e=(1,-1,-1,1)^\top$, so $SE=6/(5\sqrt{10})$.

In the conditional benchmark, $M_VD_xb=(-3,3,3,-3)^\top/(2\sqrt{10})$. The numerator and denominator weight columns are respectively even and odd under site reversal, so their projected Gaussian coordinates are independent. The same parity holds in the location model. Their ratio at $q=1$ is Cauchy, with scale given by the square root of the variance ratio.

In the location benchmark the numerator variance is $A=1+1.5r+r^2+0.5r^3$ and the denominator variance is $B=1+0.6r-0.8r^2-0.8r^3$. On the admissible interval $B>0$ and $A'B-AB'=0.9+3.6r+5.7r^2+3r^3+0.4r^4>0$. In the conditional benchmark the numerator variance $2.5+2.5r+2r^2+2r^3$ increases and the denominator variance $(9/40)(4-6r+4r^2-2r^3)$ decreases. Both rejection suprema are therefore at $r=1/2$. At this endpoint, the variances are $(33/16,1)$ for location and $(9/2,63/160)$ for conditional OLS. For each Cauchy scale $s$, the tail at $k$ is $2\operatorname{atan2}(s,k)/\pi$ and the cutoff is $s\cot(\pi\alpha/2)$.

This yields a conditional p-value of $0.5784795418841093$ and the combined interval $[-15.300205077361184,17.300205077361184]$. With the historical order $D_xM_VS_xb$, the denominator column is instead $(1,-1/2,1/2,-1)^\top/\sqrt{10}$ up to sign. Its variance is $3/20$ at $c_{\min}$, giving a conditional tail of $0.714517989950016$ at that decay and the same observed t-statistic. The fixture checks this fixed-decay contrast directly. [Conditional SCPC definition](../scpc-conditional-ols.md#two-gaussian-benchmark-families)

The historical order can also give too small a cutoff. Keep the four sites, $q=1$, $\bar\rho_{\max}=17/48$, and the intercept, and take $x=(-1,-3/2,3/2,1)^\top$. Site reversal still makes the numerator and denominator independent. For both orders the numerator variance increases and the denominator variance decreases over $0\le r\le1/2$. At $r=1/2$ the numerator variance is $27/8$, with denominator variance $7/640$ in the paper model and $1/16$ under the historical order. The nominal 5% cutoffs are about $223.1998$ and $93.3712$, both above the location cutoff. Evaluating the historical cutoff under the paper null at this decay gives rejection probability $0.11838546$.

## SCPC component selection

Take three equilateral sites with unit edge lengths, an intercept-only fit to $y=(1,2,3)^\top$, and $\bar\rho_{\max}=1/4$. For average correlation $r$, the constant direction of the Gaussian covariance has variance $1+2r$. Every centered unit vector has variance $1-r$. These coordinates are jointly independent, and the two benchmark families coincide for this design. The statistic has law $\sqrt{(1+2r)/(1-r)}\,t_q$, with maximal tails at $r=1/4$. Compound symmetry makes each candidate law independent of the choice inside the tied centered eigenspace.

At 95%, the two candidate critical values are about $17.9693$ and $6.08487$. Multiplying by $2E[\sqrt{\chi^2_q/q}]/\sqrt n$ gives i.i.d. expected lengths of about $16.5554$ and $6.22681$. Hence $q=2$. The full centered space then gives coefficient $2$, $SE=1/\sqrt3$, and p-value $1-\sqrt3/2$ for the zero null. [Component selection definition](../scpc-conditional-ols.md#recommended-choice-of-component-count)

## SPUR half-life continuous inversion

Take locations $0,1,2$, $q=2$ and $y=(1,-6,0)^\top$. The centered Lévy–Brownian covariance has eigenvalues $1,1/3$ with unit directions $(1,0,-1)^\top/\sqrt2$ and $(1,-2,1)^\top/\sqrt6$. For $r=2^{-1/h}$, the projected exponential covariance is diagonal, with entries $1-r^2$ and $(1-r)(3-r)/3$ before the common $1/(2c)$ factor. Its scale-free variance ratio is

$$
\kappa(h)=\frac{3(1+r)}{3-r}.
$$

The first squared directional share is $u=z_1^2/(z_1^2+z_2^2)=3/172$. For any $\kappa$, the angular density at that direction is proportional to $\sqrt{\kappa}/(\kappa+(1-\kappa)u)$. Integrating over $h\in[0,2]$ with uniform weight yields the normalized mixture density $0.7898748172656417$ with the dimension-only constant omitted. Dividing by the null density at $h_0=2$ gives $L_2=1.169105060364498$.

For $h_0\ge 2$, each alternative has $\kappa(h)\le\kappa(h_0)$, so the ratio is strictly decreasing in $u$. The Gaussian angular tail therefore reduces to a Cauchy probability:

$$
p(h_0)=\frac2\pi\arctan\!\left(\frac{\sqrt3}{13\sqrt{\kappa(h_0)}}\right).
$$

The unique solution of $p(h_0)=0.05$ is $14.814779473289729$. It follows by solving $\kappa_*=[\sqrt3/(13\tan(0.05\pi/2))]^2$, $r_*=3(\kappa_*-1)/(\kappa_*+3)$, and $h_*=-\log2/\log r_*$.

For the remaining nulls $h_0\in[0,2]$, write $A(\kappa,u)=\kappa+(1-\kappa)u$, $\kappa_0=\kappa(h_0)$, and $T(u)=E_h[\sqrt{\kappa/\kappa_0}A(\kappa_0,u)/A(\kappa,u)]$. The derivative satisfies

$$
\frac{d}{du}\left[A(\kappa_0,u)^2T'(u)\right]=
2A(\kappa_0,u)\,E_h\left[\frac{\sqrt{\kappa/\kappa_0}(\kappa-\kappa_0)^2}{A(\kappa,u)^3}\right]>0.
$$

If $T'$ at the observed share is nonpositive, all smaller shares give a larger test statistic, and their probability is at least $p(2)>0.05$. If $T'$ is positive, all larger shares give a larger statistic, and their probability is at least $1-(2/\pi)\arctan(\sqrt3/13)>0.9$. So all these nulls are retained, and the confidence-set closure is $[0,h_*]$ in distance units. This also certifies the unbounded variant $y=(1,-5,0)^\top$. Its ratio of observed directions is $\sqrt3/11$, giving $p(\infty)=0.057715876752608954$, and the same argument retains all smaller half-lives. [Half-life definition](../spur-half-life.md#confidence-set-inversion)
