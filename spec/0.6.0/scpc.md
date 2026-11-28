# SCPC: statistical definition and supported numerical procedure

This chapter is normative for methods **0.6.0**, together with the [scope and geometry contract](README.md) and [numerical profile](numerics.md). It defines the adopted target; the baseline implementation remains methods 0.5.0 until the [conformance checks](verification.md) pass. No historical specification is needed to determine the formulas below.

Provenance labels are **[P]** published mathematical definition/result, **[R]** published recommendation, **[A]** author-software convention, **[D]** package decision, and **[I]** algebraic consequence. A numerical choice can implement a paper definition without itself being a theorem. Source identities and author-software pins are recorded in [sources.md](sources.md).

## 1. Observations, sites, and fitted-model contract

**[D]** Let $N$ be the number of fitted observations and $m\ge2$ the number of declared spatial sites. All observation-level arrays use the retained fitted order. Let $A\in\{0,1\}^{N\times m}$ have exactly one unit entry in each row, mapping each fitted observation to its declared site. Every site has at least one observation and one coordinate vector. Without explicit site clustering, $m=N$ and $A=I_N$. Different declared sites may have identical coordinates; this does not merge their fitted observations or redefine their regression weights.

The public supported fits are those in [README.md](README.md): ordinary unweighted OLS and the specified ordinary unweighted 2SLS fits, with the retained fitted nuisance structure. This chapter also covers the mean as an intercept-only OLS fit. A class name, a residual vector, or a successful numerical rank check does not by itself establish the estimator contract. Penalized fits and estimator-changing subclasses must not be interpreted as OLS. Unsupported weights, offsets, and estimators receive the specified errors.

Write $X\in\mathbb R^{N\times p}$ for the full retained structural design and $y\in\mathbb R^N$ for the fitted outcome. $X$ includes all fitted, nonabsorbed nuisance regressors, not only coefficients requested for display. When an intercept is fitted explicitly, it is a column of $X$. An explicit no-intercept OLS fit is permissible for SCPC; it must retain its actual design.

Let $F$ be the orthogonal residual-maker for the complete absorbed fixed-effect space. Set $F=I_N$ if no effects are absorbed. Define $X_w=FX$ and $y_w=Fy$. Host within-transformed matrices and residualizers must represent this same fitted space and sample. The coefficient selector $e_j\in\mathbb R^p$ has a unit entry at the requested coefficient; it is not a regression disturbance.

### 1.1 Ordinary OLS

**[P/I]** Define the bread, influence direction, and full residual-maker by

$$
B=(X_w^\top X_w)^{-1},\qquad
a_j=X_wBe_j,\qquad
M=F-X_wBX_w^\top.
$$

Then $\widehat\beta_j=a_j^\top y$ and $\widehat e=My$. These are full-design expressions for the Frisch–Waugh–Lovell coefficient and its fitted residuals. If $x_j$ is the target regressor residualized against all other fitted regressors and effects, then $a_j=x_j/(x_j^\top x_j)$. The stored fitted estimates and residuals must agree with these estimator definitions within the host-validation tolerances.

The bread is the inverse **total** cross product. It contains no residual variance, degrees-of-freedom adjustment, or additional factor of $N$ or $m$. Per-observation scores may be supplied as $X_w\odot\widehat e$; multiplying their coefficient columns by $B$ gives the influence contributions $a_j\odot\widehat e$. See [Müller and Watson (2023), §§2.1–2.2](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf).

### 1.2 Ordinary 2SLS

**[P/I/D]** Let $Z$ contain the full instrument design, including fitted exogenous regressors and the intercept when applicable. Set $Z_w=FZ$, and let $P_Z$ be the orthogonal projector onto the column space of $Z_w$. Redundant instrument columns do not change $P_Z$; the projected structural design must identify the retained coefficients.

$$
B=(X_w^\top P_ZX_w)^{-1},\qquad
a_j=P_ZX_wBe_j,\qquad
M=F-X_wBX_w^\top P_Z.
$$

Again $\widehat\beta_j=a_j^\top y$ and $\widehat e=My$. The score design is $P_ZX_w$, not the raw endogenous regressor. Its coefficients, columns, instruments, fitted transformations, and sample must be aligned before calibration. The inverse notation defines the mathematics; implementations must use stable solves/factorizations under [numerics.md](numerics.md).

The residual-maker is generally oblique:

$$
M^\top=F-P_ZX_wBX_w^\top\ne M.
$$

Strong identification and the relevant score-limit assumptions are required for the ordinary IV interpretation. Conditional finite-sample Gaussian reasoning that also conditions on endogenous structural regressors is heuristic under endogeneity. The conditional cutoff does not repair weak identification. [Müller and Watson (2023), §5.3](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf) discusses the scalar Anderson–Rubin alternative and the limitation of conditioning on endogenous regressors. The package support scope is in [README.md](README.md).

### 1.3 Site influences and sample sizes

**[P/I]** For each requested coefficient, form the site influence vector

$$
h_j=A^\top(a_j\odot\widehat e)\in\mathbb R^m.
$$

This is a sum of observation influences within a declared site, not an average. Clustering does not refit the model, replace observation weights, or substitute site count for observation count in the fitted cross products. The coefficient calculation remains at the observation level; the spatial basis and its normalization use the $m$ declared sites.

## 2. Spatial benchmark and exact endpoint

**[P]** For the physical site distance matrix $D$ and positive sample diameter $\Delta=\max_{g,h}D_{gh}$ defined in [README.md](README.md#52-metrics-and-units), define normalized site distances $d=D/\Delta\in[0,1]$ and

$$
K(c)_{gh}=\exp(-c d_{gh}),\qquad
\bar\rho(c)=\frac{1}{m(m-1)}\sum_{g\ne h}K(c)_{gh}.
$$

Internal decay rates $c$ ($c_0$, $c_{\rm end}$, and mesh points) multiply normalized distances $d$; public physical-rate views convert via $c_{\rm physical}=c/\Delta$ as stated in §8. The correlation bound is $\mathrm{avc}$, and $c_0$ solves $\bar\rho(c_0)=\mathrm{avc}$. The default (`0.03`) and calibration rules are fixed in [numerics.md](numerics.md). The public domain is the mathematical attainability interval below; there is no separate historical `0.001`/`0.99` cutoff. The exponential benchmark and calibration by average pairwise correlation follow [Müller and Watson (2022), §2 and Remark 2.3](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf).

**[I/D]** Define the co-location matrix $J$ and off-diagonal co-location fraction $f$ by

$$
J_{gh}=\mathbf1\{d_{gh}=0\},\qquad
f=\frac{1}{m(m-1)}\sum_{g\ne h}J_{gh}.
$$

For a geometry containing a positive distance, $\bar\rho(c)$ decreases strictly from $1$ to $f$. A finite positive $c_0$ requires $f<\mathrm{avc}<1$. Wholly co-located data do not meet this condition. Equality at $f$, and unattainable or numerically unresolved decay targets, must produce the specified calibration error.

The exact infinite-decay limit is

$$
K(\infty)=J.
$$

No diagonal nugget or additional IID covariance is introduced. With distinct site coordinates $J=I_m$; with duplicates it is block-valued. This limit follows from the selected kernel. Support for distinct declared sites with identical coordinates is a package decision.

For any weight matrix $W\in\mathbb R^{m\times r}$, compute its endpoint covariance as

$$
W^\top JW
=\sum_{\ell}
\left(\sum_{g:\,s_g=s_\ell}W_{g,:}\right)^\top
\left(\sum_{g:\,s_g=s_\ell}W_{g,:}\right),
$$

where the outer sum ranges over distinct coordinates. Replacing this expression by $W^\top W$ when coordinates coincide changes the benchmark.

**[D]** The practical finite mesh $\mathcal C$, consisting of the single uniform log-rate grid from $c_0$ to $c_{\rm end}$ and the exact limit $K(\infty)$, is defined by [numerics.md](numerics.md). The same full-site geometry and mesh must serve exact and Nyström calculations, every candidate $q$, every coefficient, every reported significance level, and every p-value. Full-site decay limits and the mesh are computed once per operation, outside component expansion.

## 3. Working spatial basis and statistic

### 3.1 Exact basis

**[P/I]** Let

$$
b_0=\mathbf1_m/\sqrt m,\qquad
C_m=I_m-b_0b_0^\top,\qquad
H_0=C_mK(c_0)C_m.
$$

For an eligible $q$, $R_q\in\mathbb R^{m\times q}$ contains the leading positive eigenvectors of $H_0$, normalized so that

$$
R_q^\top R_q=I_q,\qquad R_q^\top\mathbf1_m=0.
$$

The paper's eigenvectors have norm $\sqrt m$; this chapter uses unit-norm columns throughout. The normalized constant is stored separately in $W_{\rm loc,q}=[b_0,R_q]$. It is not one of the $q$ variance directions. See [Müller and Watson (2022), equations (3)–(6)](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf).

**[D]** Numerical positive rank and eigenvalue blocks follow [numerics.md](numerics.md). Candidate prefixes must end at complete positive blocks. A cutoff inside a numerical tie must be extended to a complete block or rejected as specified; it must not be resolved by an arbitrary provider eigenvector or a randomized partial block. Rotations within a fully retained subspace leave the SCPC statistic and covariance-calibrated probabilities unchanged. Exact individual eigenvector entries are not the invariant output.

### 3.2 Standard error and observed statistic

**[P/I]** For coefficient $j$ and scalar null $\beta_{j,0}$, use

$$
SE_j=\sqrt{\frac{m}{q}}\,\|R_q^\top h_j\|_2,
\qquad
T_j(\beta_{j,0})=\frac{\widehat\beta_j-\beta_{j,0}}{SE_j}.
$$

Because $R_q^\top\mathbf1_m=0$, centering $h_j$ before its projection is an algebraically equivalent numerical stabilization. There is no residual-degrees-of-freedom multiplier. Internal score arrays scaled by $m$ are permissible only if the compensating standard-error normalization gives this formula.

Reported ordinary coefficient p-values use $\beta_{j,0}=0$. If an adapter exposes a separate scalar-hypothesis calculation, it uses the same $SE_j$, working basis, and calibration family with its stated null; this definition does not require an additional public testing interface. Every interval is marginal; a collection of coefficient intervals is not a simultaneous confidence region.

## 4. The two calibration experiments

### 4.1 Location family

**[P]** The location experiment is a Gaussian site-score vector with covariance $K(c)$, projected through

$$
W_{\rm loc,q}=[b_0,R_q],\qquad
\Omega_{\rm loc,q}(c)=W_{\rm loc,q}^\top K(c)W_{\rm loc,q}.
$$

This describes the calibration experiment. Applying its cutoff to a fitted-regression score statistic has the paper's applicable asymptotic assumptions; it is not a finite-sample assertion for every fitted regressor or arbitrary heterogeneous score process. Unconditional SCPC uses this family alone.

### 4.2 Coefficient-specific conditional family

**[P/I]** Define

$$
\Xi_j=\operatorname{diag}(a_j)A\in\mathbb R^{N\times m}.
$$

For each site with a nonzero column of $\Xi_j$, normalize that column to unit Euclidean length to obtain $S_j\in\mathbb R^{N\times m}$. Zero columns use the explicit singleton/cluster convention below; in particular the singleton convention is $\operatorname{sgn}_+(0)=+1$. The conditional benchmark is

$$
\eta\mid(X,Z,F,A,s)\sim N(0,K(c)),\qquad
e=S_j\eta,\qquad \widehat e=M S_j\eta.
$$

This is the specified conditional site-direction model. In particular it uses a rank-one disturbance direction within each declared site with nonzero score direction; accepting a cluster label does not assert finite-sample protection for every within-site covariance matrix.

The weight columns are

$$
W_{{\rm cond},j,q}
=\left[S_j^\top\Xi_j b_0,
\;S_j^\top M^\top\Xi_j R_q\right],
\qquad
\Omega_{{\rm cond},j,q}(c)
=W_{{\rm cond},j,q}^\top K(c)W_{{\rm cond},j,q}.
$$

The first projected coordinate equals $(\widehat\beta_j-\beta_{j,0})/\sqrt m$ under the benchmark null; the remaining coordinates generate the variance projection. The formulas transpose the observation-level numerator and denominator maps. See [Müller and Watson (2023), §§2.3 and 5.1](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf).

For one observation per site, write $D_a=\operatorname{diag}(a_j)$ and $S_a=\operatorname{diag}(\operatorname{sgn}_+(a_j))$. Then

$$
W_{{\rm cond},j,q}=[S_aD_ab_0,\;S_aM^\top D_aR_q].
$$

For OLS, $M^\top=M$; for IV the transpose must be retained. Multiplication by the direction and fitted residualization do not commute.

**[D/A] Zero directions.** The adopted singleton convention is $\operatorname{sgn}_+(0)=+1$, matching the stated author MATLAB convention. For a site containing multiple observations whose entire coefficient direction is zero, its column of $S_j$ is the zero vector. This explicitly sets that site's disturbance direction to zero in this particular conditional benchmark. It is not a unique limiting direction as the score tends to zero and is not a claim that the observed disturbances there vanish. The location family remains part of conditional SCPC. These conventions must not vary silently across hosts.

### 4.3 Fixed effects and IV operator order

**[I]** To obtain the conditional denominator weights without materializing an $N\times N$ matrix:

1. Form the observation-level block $\Xi_jR_q$.
2. Apply the full fitted fixed-effect residualizer $F$.
3. Apply $I_N-X_wBX_w^\top$ for OLS, or $I_N-P_ZX_wBX_w^\top$ for IV.
4. Multiply by $S_j^\top$, aggregating to declared sites.

The products in steps 2–3 equal $M^\top$ because the within designs lie in the range of $F$. Applying only the within-slope projection omits the fixed effects. Aggregating to sites before these observation-level operations generally changes the result.

A common nonzero rescaling of $a_j$, applied consistently in the conditional construction, does not change the Gaussian ratio probabilities. Such normalization must not alter the physical units of the separately reported coefficient or standard error.

When $P_ZX_w=X_w$, the IV formula reduces to OLS. Conditional weights, probabilities, and intervals must then agree for identical fitted influence contributions and geometry, within numerical tolerances. This is an algebraic conformance check, not a requirement that every host estimator accept every identity-instrument fit.

## 5. Common family, critical values, and p-values

**[P/I]** For a projected covariance $\Omega$, let $G\sim N(0,\Omega)$ and define the strict rejection probability in ordinary $t$ units by

$$
r_\Omega(k)=\Pr\!\left(G_0^2-\frac{k^2}{q}\sum_{\ell=1}^{q}G_\ell^2>0\right),
\qquad k\ge0.
$$

The quadratic-form event remains defined for singular positive semidefinite $\Omega$; implementations must not divide by a simulated zero denominator. Covariance validation and degenerate boundary evaluation follow [numerics.md](numerics.md).

For unconditional inference define

$$
\mathcal F_{j,q}=\{\Omega_{\rm loc,q}(c):c\in\mathcal C\}.
$$

For conditional inference define

$$
\mathcal F_{j,q}
=\{\Omega_{\rm loc,q}(c),\Omega_{{\rm cond},j,q}(c):c\in\mathcal C\}.
$$

Let $R_{j,q}(k)=\max_{\Omega\in\mathcal F_{j,q}}r_\Omega(k)$. The practical critical-value target is

$$
k_{j,q,\alpha}=\inf\{k\ge0:R_{j,q}(k)\le\alpha\}.
$$

The solver returns the conservative bracket endpoint under the critical-value and probability tolerances in [numerics.md](numerics.md). The same family and mesh are used for the observed two-sided p-value. Define its observed tail inclusively:

$$
p_{j,q}(t)=\max_{\Omega\in\mathcal F_{j,q}}
\Pr\!\left(G_0^2-\frac{t^2}{q}\sum_{\ell=1}^{q}G_\ell^2\ge0\right),
\qquad t=|T_j|.
$$

**[D/I]** For continuous ratio laws this equals $R_{j,q}(|T_j|)$. At a genuine atom, the equality probability must be included; silently using the strict tail can report a zero p-value at the constant value of a degenerate ratio. The cutoff still controls the strict rejection event. Closed intervals retain equality at their endpoint, and the operational interval decision uses $|T_j|>k_{j,q,\alpha}$. For continuous laws, usual p-value/interval comparisons can differ only at the exact nominal boundary or within the stated numerical tolerance. Neither rounded printed p-values nor a strict-tail substitution at an atom overrides that convention.

Genuine atom cases may be evaluated algebraically from the projected Gaussian experiment. If a required degenerate boundary cannot be resolved within the supported numerical precision, report a clear calibration-degeneracy error instead of inventing finite significance or elaborate unverified fallback behavior. The precision and equality-detection rules are in [numerics.md](numerics.md).

**[I/D] No finite cutoff.** Write $\Omega_{DD}$ for the covariance block of the $q$ denominator coordinates. If a required valid covariance has $\Omega_{00}>0$ and $\Omega_{DD}=0$, positive semidefiniteness also makes its numerator-denominator cross-covariances zero. Then $r_\Omega(k)=1$ for every finite $k$, and the extended-real critical-value target is $+\infty$; the inclusive tail is likewise one for every finite observed $t$. This is algebraic degeneracy, not failure to find a sufficiently large numerical bracket.

The supported profile reports a clear calibration-degeneracy error for the affected coefficient, identifying the family and covariance point and explaining that the benchmark leaves no random estimated-variance direction. It does not add a special successful infinite-interval path. For conditional SCPC, keep the location-selected $q$: do not reselect it, drop the offending covariance, or repeatedly expand a root bracket to conceal this condition. Co-location-group cancellation can cause it even when the unprojected conditional weight norm passes the projection check, and can hold at finite rates as well as at $K(\infty)$. Check the covariance after projection onto the actual kernel support. An exactly zero numerator and denominator law is separately undefined for the ratio and also receives a degeneracy error. A numerically unresolved small covariance is a precision error; a tolerance comparison alone is not proof of an exactly infinite mathematical cutoff.

The marginal confidence interval is

$$
[\widehat\beta_j-k_{j,q,\alpha}SE_j,
\;\widehat\beta_j+k_{j,q,\alpha}SE_j].
$$

**[D] No inherited Student floor.** Initialize bracketing from a simple positive starting guess or previously evaluated cutoff, and never impose a Student-$t_q$ lower floor; no Student-$t$ quantile implementation is required. Bracketing must be allowed downward when required. With distinct site coordinates the IID location benchmark is already $K(\infty)$. With duplicate coordinates it is not an additional member of the adopted family. Neither critical values nor p-values may silently add IID, a nugget, or a historical compatibility floor.

Every point in the single finite-rate mesh and exact endpoint required by [numerics.md](numerics.md) enters both the cutoff and p-value calculations. Conditional intervals use the larger cutoff required by the two families; conditional p-values use their common worst-case tail definition. Reporting a location-only p-value beside a conditional interval is not this method.

**[P/D] Continuous target versus shipped approximation.** The paper target replaces the finite maximum by a supremum over all $c\ge c_0$, with the stated limit. Methods 0.6.0 adopts the single bounded log-mesh procedure and accuracy status in [numerics.md](numerics.md). Extending the mesh across short-scale correlations to $36/d_{\min}$ addresses the demonstrated omitted-tail case without certifying the entire continuum. “Exact SCPC” means full spatial geometry and eigenspace calculation, not exact continuous nuisance maximization. The returned numerical status must preserve this distinction.

## 6. Choosing the component count

**[P/R]** Choose $q$ using the location family at size $0.05$, independently of the observed outcomes and coefficient-specific conditional families. For each eligible candidate define

$$
L(q)=k_{{\rm loc},q,0.05}
\frac{\Gamma((q+1)/2)}{\sqrt q\,\Gamma(q/2)}.
$$

For integer $q\ge1$, $r_q=\Gamma((q+1)/2)/\Gamma(q/2)$ satisfies the exact two-step recurrence $r_1=1/\sqrt{\pi}$, $r_2=\sqrt{\pi}/2$, $r_{q+2}=\frac{q+1}{q}r_q$, requiring no general log-gamma routine. This objective is proportional to expected interval length under a separate Gaussian IID site-score experiment, because the columns of $R_q$ are orthonormal. That reference experiment is an efficiency criterion; it does **not** add IID to the calibration family when coordinates coincide. Minimize $L(q)$ over the evaluated eligible candidates, breaking exact objective ties toward smaller $q$.

Using the 5% selection level and retaining the same $q$ for all reported confidence levels and p-values follows [Müller and Watson (2022), equation (6) and Remark 2.6](https://www.princeton.edu/~mwatson/papers/ECTA19465.pdf). Conditional calibration is performed after this choice. It must not reselect $q$ separately for whichever coefficient or null happens to produce a preferred result.

**[A/D] Finite candidate search.** The finite candidate-prefix rule, expansion by doubling, and complete-block stopping rule are governed normatively by [numerics.md](numerics.md#23-scpc-candidate-expansion-and-landmark-sampler) (starting at $\min(\mathrm{limit},20)$ and doubling $k\leftarrow\min(\mathrm{limit},2k)$). A single Nyström sample does not disable expansion.

Stopping at an interior winner is an inherited computational heuristic, not proof of the global minimum over every positive dimension. A nonoptimal $q$ affects efficiency; calibration must still use the actual chosen weights and common covariance family. This specification does not require an exhaustive search merely to preserve an unwarranted optimality claim.

## 7. Nyström working basis

### 7.1 Target and sampling

**[R/D]** Landmark extension and pooling are motivated by [Müller and Watson (2023), §5.2, equation (10)](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf). The package's particular weighting and pooling below define its approximate working basis; they are not a verbatim implementation of that equation or the author software.

Defaults are 1,000 landmark sites, 20 samples, and landmark seed 1. Use $\ell=\min(1000,m)$, or the validated requested landmark count. Samples are drawn without replacement within each subset from canonical site order, using the versioned sampler and seed contract in [numerics.md](numerics.md). Subset choices depend on geometry and the supplied seed, never outcomes, fitted residuals, coefficient significance, or target-specific leverage. Reuse the same subsets when expanding the working dimension.

Canonical geographic order uses latitude then longitude after the common coordinate canonicalization; Euclidean order is lexicographic. Numerical translation/common scaling may condition Euclidean coordinates without changing distance ratios. Restore weight rows to reported site order. If $\ell=m$, use a single complete sample and return the exact working subspace, allowing an exact-algorithm shortcut.

The complete full-site geometry determines $c_0$, the finite mesh limits, and $\mathcal C$. Landmark-pair averages do not replace those calculations. Every selected landmark subset must yield positive rank and finite nonzero embedded weights; if any landmark subset has zero positive rank or no usable embedded norm, fail clearly with an actionable error instructing the caller to increase the landmark count or use exact geometry—do not skip unusable subsets, silently continue with survivors, resample subsets automatically, or alter the global correlation target.

### 7.2 Embedded and pooled operator

**[D/I]** For subset $S_r$ of size $\ell$, obtain complete leading positive blocks $U_r$ of its centered covariance $H_r=C_\ell K_{S_r,S_r}(c_0)C_\ell$ at the common $c_0$, where $C_\ell=I_\ell-\mathbf1_\ell\mathbf1_\ell^\top/\ell$. Apply the unified spatial rank and block-gap rule in [numerics.md](numerics.md#21-rank-and-eigenvalue-blocks) using operator size $\ell$ and this local matrix's own $\lambda_{\max}$; the pooled operator below uses operator size $m$ and its own $\lambda_{\max}$, so the two scales are not interchanged. With $A_r=U_r^\top H_r U_r$ (which equals the diagonal eigenvalue matrix $\Lambda_r$ in an eigenbasis, using actual retained eigenvalues without block-mean averaging), define

$$
E_r=C_mK_{:,S_r}(c_0)U_rA_r^{-1/2},
\qquad
\widetilde E_r=E_r/\|E_r\|_F.
$$

Under any orthogonal rotation $U_r\mapsto U_rO$ of complete retained blocks, $A_r^{-1/2}\mapsto O^\top A_r^{-1/2}O$ so $E_r\mapsto E_rO$ and $E_rE_r^\top$ is invariant without artificial eigenvalue averaging. Use one Frobenius scale for the whole embedded sample; per-column normalization is not substituted because it breaks rotation invariance within a tied local block. Any sample with zero or nonfinite embedded norm produces the clear failure in §7.1.

Concatenate the retained sample blocks,

$$
E=[\widetilde E_1,\ldots,\widetilde E_J],\qquad H_{\rm pool}=EE^\top.
$$

The working $R_q$ is the leading complete positive eigenspace of $H_{\rm pool}$, centered and orthonormalized as required by [numerics.md](numerics.md). This defines the target even when an implementation uses the smaller Gram matrix $E^\top E$ rather than materializing $H_{\rm pool}$. Forming a Gram matrix squares its condition number; final centering/orthogonality and retained-rank checks must be performed on the recovered site directions.

Dimension expansion may change the local retained blocks and hence the pooled operator. Recompute the candidate objective for the actual resulting basis. The admissible working dimension is at most $\ell-1$, as well as at most $m-1$ and the available positive pooled rank. If a required complete boundary block cannot be represented within the requested landmark limit, report that limit and recommend more landmarks or exact geometry instead of splitting the block.

### 7.3 Calibration, approximation status, and cost

**[D/I]** Project every full-site $K(c)$, including $J$, onto the actual selected weights. Random sampling of covariance pairs is not part of this profile. Approximate principal components do not by themselves invalidate Gaussian calibration of that statistic over the evaluated family. Finite landmark counts do not guarantee equivalence to exact principal components, identical $q$, or identical efficiency.

Full-pair geometry/calibration still requires quadratic pair work in $m$. Streaming decay calculations in $O(m)$ auxiliary memory does not make the complete operation linear-time or describe all of its allocations. Landmark decompositions, the embedded $m$-by-total-directions array, pooled Gram work, and covariance projections have their own costs. Reuse geometry-only limits across $q$ expansion; no persistent global cache is required.

The current automatic policy selects Nyström at 4,500 sites and exact below that threshold. This is a package memory/runtime tradeoff, not a theorem or a universal speed crossover. A forced exact or Nyström request remains explicit. A result records the realized method, sample count and landmark count, seed/sampler identity, and approximation status. Sensitivity to landmark settings can be examined without changing the stated statistical family.

## 8. Degeneracy, units, and result contract

**[D]** Fail clearly for an unsupported estimator, unidentifiable coefficient, invalid site mapping, unattainable decay target, no complete positive spatial direction, or a numerically degenerate observed variance projection. A zero estimated standard error gives an explicit degenerate-inference outcome with no ordinary $t$, p-value, or confidence interval, or a clear error. It is not converted into a zero-width interval or manufactured significance. Rank and collapse tests use the scale-aware criteria in [numerics.md](numerics.md).

Singular positive semidefinite projected covariances are not automatically errors; their Gaussian quadratic-form events are evaluated with the specified boundary convention. The proven no-finite-cutoff case in section 5 receives its specific calibration-degeneracy error before root finding. A materially indefinite covariance, an unbracketed otherwise finite cutoff, or unmet required numerical checks must produce the profile's explicit error/status, rather than an apparently successful finite interval. Reaching a computational cap does not prove that a mathematical finite cutoff does not exist.

Decay rates must have recorded units. Internally, covariance calculations normalize distances by the positive sample diameter $\Delta$ ($d=D/\Delta$), and convert public decay rates via $c_{\mathrm{physical}}=c_{\mathrm{normalized}}/\Delta$. Geographic SCPC decay rates are reported in inverse meters (`1/m`, using the great-circle arc metric on radius 6,371,008.8 meters), consistent with geographic half-life (meters) and ISO radius (meters); Euclidean SCPC rates are reported in inverse input coordinate units. In either case $K(c)$ and inference are invariant to a common positive change of distance units when rates are transformed consistently.

The primary result contains coefficient identity and estimate, $SE_j$, zero-null $T_j$, the combined two-sided p-value, the 95% marginal interval, and its effective coefficient-specific critical value. Optional critical values correspond to two-sided sizes 32%, 10%, 5%, and 1% (coverages 68%, 90%, 95%, and 99%) using the same selected $q$. A separately exposed location-only cutoff must be labeled as such and must not replace the coefficient-specific cutoff in conditional intervals.

Result metadata identifies observation/site counts and order, membership, coordinate units and scale, $\mathrm{avc}$, co-location fraction, $c_0$, $c_{\rm end}$ and the evaluated mesh status, $q$, complete-block/tie status, inference mode, approximation settings, numerical tolerances/status, and methods version. Finite calibration over the selected mesh must not be labeled a certified continuous-family guarantee.

The required identities and regression examples are specified in [verification.md](verification.md). They include duplicate-site p-value/interval consistency after removing the Student floor, genuine-atom inclusive tails or the declared degeneracy error, zero-standard-error handling, conditional-tail calibration, OLS/IV equivalence, unit changes, site aggregation, and exact/full-landmark subspace agreement. Atom checks must distinguish a strict cutoff rejection probability from the inclusive observed p-value rather than asserting that those two quantities are always identical.
