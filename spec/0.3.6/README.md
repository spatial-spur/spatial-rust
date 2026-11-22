# SCPC and SPUR methods specification 0.3.6

This record extends [0.3.5](../0.3.5/README.md).

## SPUR residual response

For the evaluated full-rank design $X$, apply the residual projection $M_X$ to the response. This gives $M_X(y + X\beta) = M_Xy$ in exact arithmetic. An exact constant column permits subtracting a common outcome offset before rescaling and residualizing. This preserves small representable differences near a large offset.

The numerical constant-span check bounds the squared residual norm of the unit constant by machine epsilon. That check controls the size of the constant covariance term under projection. For a design without an exact constant column, project the scaled original outcome onto the residual space. Subtracting $y_1\mathbf{1}$ in that case can introduce a term proportional to $y_1 M_X\mathbf{1}$.

The regression test uses 60 equally spaced sites, $x_i=1+10^{-8}s_i$, and $y_i=u_i+\beta x_i$. The residual statistic agrees for $\beta=0$ and $\beta=10^9$ to floating-point tolerance. [Linked review and reproducer](https://chatgpt.com/s/t_6b02f94aa44481918787f9df5b83c974)

## Native methods compatibility

The host packages check the native API and require methods spec 0.3.6 or newer before dispatch. Auto selection uses the host fallback when a build lacks the required method version. Explicit Rust selection reports the required build. The wire contracts use SCPC API 7 and SPUR API 4.
