# Conditional IV SCPC projection

Methods specification 0.2.0. This extends the required OLS weight order from version 0.1.0 to the IV residual maker in Müller and Watson (2023), §5.3. Model extraction, coefficient naming, and observation alignment are package responsibilities.

Let $X$ contain the full structural regressor design and $Z$ the full instrument design, including the fitted exogenous regressors and the intercept when fitted. With $P_Z$ the orthogonal projector onto the instrument space, set

$$
H_{IV}=X(X^\top P_ZX)^{-1}X^\top P_Z,
\qquad M_{IV}=I-H_{IV}.
$$

The regular case requires $X^\top P_ZX$ to be nonsingular. The residual maker $M_{IV}$ is generally oblique: $M_{IV}^\top\ne M_{IV}$. An orthogonal OLS projection cannot be substituted for it.

For coefficient index $j$, the relevant score direction is proportional to $P_Z X(X^\top P_ZX)^{-\top} e_j$, with the corresponding within-design expression when fixed effects are absorbed. Here $e_j$ is the coefficient selector, not the error vector. This direction is generally different from the structural endogenous regressor. A common nonzero scale cancels from the statistic and calibration.

Use the paper's cluster maps $\Xi$ and $\Xi_s$, constructed from the coefficient's score direction, with $\Xi_s$ normalizing each nonzero cluster column. Under the conditional model $e=\Xi_s a$, fitted errors satisfy $\widehat e^{IV}=M_{IV}\Xi_s a$. For cluster-level spatial basis $R$, the denominator projections are

$$
R^\top\Xi^\top\widehat e^{IV}
=R^\top\Xi^\top M_{IV}\Xi_s a.
$$

The corresponding weight **columns**, used to form $W^\top K(c)W$, must therefore be

$$
W_{\mathrm{den}}=\Xi_s^\top M_{IV}^\top\Xi R.
$$

For one observation per cluster, these reduce to $S_x M_{IV}^\top D_xR$. Multiply by the ordinary regressor/score direction before projection, apply the **transposed** IV residual maker, and multiply by the normalized/sign direction after projection. The numerator column remains $\Xi_s^\top\Xi b_0$. The conditional interval uses the same combined family cutoff convention as OLS.

If $M_{IV}=I-UV^\top$, a generic operation $A\mapsto A-U(V^\top A)$ applies $M_{IV}^\top$ by taking its factor inputs as $V,U$. This factor representation is an implementation choice. The transposed weight columns are the mathematical requirement for every backend.

For absorbed fixed effects, let $D_F$ be orthogonal projection off the full absorbed-effect space. Construct the within designs $X_w=D_F X$ and $Z_w=D_F Z$. The full residual map is $D_F-H_w$, where $H_w=X_w(X_w^\top P_{Z_w}X_w)^{-1}X_w^\top P_{Z_w}$. Its transpose is $D_F-H_w^\top$. Applying $D_F$ to the weighted basis and then applying $I-H_w^\top$ gives that transpose because the within designs are orthogonal to the absorbed-effect space. Projecting only on within slopes omits the fixed effects. All these operations remain at the observation level before cluster aggregation.

When the instrument space contains the regressor space, $P_ZX=X$, so $M_{IV}$ equals the symmetric OLS residual maker. The conditional weights, cutoff, and p-value must then agree with OLS for identical scores and spatial inputs. A host estimator may reject a perfectly explained first stage; the projection identity still provides an exact numerical oracle independent of that host policy.

The five-site identity-instrument case in the regression tests has sites $(0,1,2,4,6)$, regressor $(-1,-4,2,2,1)$, an intercept, and average-correlation bound $0.5$. Its corrected 5% cutoff is approximately $17.62678224$. The reversed IV weights previously gave $14.50861818$ and Gaussian rejection probability $0.07143375445$ at the first grid covariance. The core test checks the latter probability against independent quadrature and checks rejection control after correction.

Source: [Müller and Watson (2023), §5.3, p. 1063](https://www.princeton.edu/~mwatson/papers/Mueller_Watson_JBES_2023.pdf#page=15), together with the [version 0.1.0 source register](../0.1.0/sources.md).
