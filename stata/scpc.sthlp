{smcl}
{* *! version 0.6.0  24nov2026}{...}
{viewerjumpto "Syntax" "scpc##syntax"}{...}
{viewerjumpto "Description" "scpc##description"}{...}
{viewerjumpto "Options" "scpc##options"}{...}
{viewerjumpto "Examples" "scpc##examples"}{...}
{viewerjumpto "Stored results" "scpc##results"}{...}
{viewerjumpto "Authors" "scpc##authors"}{...}
{viewerjumpto "References" "scpc##references"}{...}
{title:Title}

{p2colset 5 14 16 2}{...}
{p2col :{bf:scpc} {hline 2}}Spatial Correlation Principal Components (SCPC) standard errors and confidence intervals{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:scpc} [{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Inference and Kernel}
{synopt :{opt k(#)}}report results for only the first {it:#} fitted coefficients{p_end}
{synopt :{opt avc(#)}}worst-case average pairwise correlation in {cmd:(0.001, 0.99)}; default is {cmd:0.03}{p_end}
{synopt :{opt cvs}}compute and display critical values at the 68%, 90%, 95%, and 99% levels{p_end}
{synopt :{opt unc:ond}}report unconditional (asymptotic normal / minimax) critical values rather than conditional SCPC critical values{p_end}
{synopt :{opt lat:long}}treat {cmd:s_1} as latitude and {cmd:s_2} as longitude (degrees) using great-circle chordal distances{p_end}
{synopt :{opt cl:uster(varname)}}group observations into spatial sites indexed by {it:varname}{p_end}

{syntab:Geometry and Approximation}
{synopt :{opt method(string)}}spatial basis method: {cmd:auto} (default), {cmd:exact}, or {cmd:nystrom} ({cmd:approx}){p_end}
{synopt :{opt exact}}alias for {cmd:method(exact)}{p_end}
{synopt :{opt approx}}alias for {cmd:method(nystrom)}{p_end}
{synopt :{opt seed(#)}}integer seed in {cmd:[0, 2^32)} for Nyström landmark selection; default is {cmd:1}{p_end}
{synopt :{opt land:marks(#)}}number of Nyström landmarks ({cmd:>= 2}); default is {cmd:1000}{p_end}
{synopt :{opt sam:ples(#)}}number of spatial averages sampled per landmark ({cmd:>= 1}); default is {cmd:20}{p_end}
{synopt :{opt tol:erance(#)}}critical-value solver tolerance ({cmd:> 0}); default is {cmd:0.001}{p_end}

{syntab:Output and Runtime}
{synopt :{opt weights}}store SCPC spatial site weights in {cmd:e(W)} and representative row indices in {cmd:e(site_rows)}{p_end}
{synopt :{opt back:end(string)}}native backend ({cmd:rust} or {cmd:auto}){p_end}
{synopt :{opt nt:hreads(#)}}thread count override for the Rust backend; {cmd:-1} uses the default runtime policy{p_end}
{synopt :{opt cap:abilities}}verify and display the loaded native plugin capabilities without modifying {cmd:e()}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:scpc} is a postestimation command implementing the Spatial Correlation Principal Components (SCPC) method of {help scpc##mw2022:Müller and Watson (2022)} and {help scpc##mw2023:Müller and Watson (2023)}.
It takes the immediately preceding unweighted {cmd:regress}, one-way {cmd:areg} (with {cmd:uncond}), or ordinary {cmd:ivregress 2sls} fit in memory, reads spatial coordinates from consecutive numeric variables {cmd:s_1}, {cmd:s_2}, ..., and reports SCPC standard errors, t-statistics, p-values, and confidence intervals while preserving the caller's regression results in {cmd:e(b)}, {cmd:e(V)}, and {cmd:e(sample)}.

{pstd}
{bf:Conditional 2SLS qualification:} For {cmd:ivregress 2sls}, valid t-test and confidence interval coverage relies on strong first-stage identification. Under weak instruments, use Anderson-Rubin inference on reduced-form projections. Under endogeneity, the finite-sample conditional Gaussian argument is heuristic, although SCPC takes the maximum over unconditional and conditional critical values.


{marker options}{...}
{title:Options}

{dlgtab:Inference and Kernel}

{phang}
{opt k(#)} restricts output and stored SCPC matrices to the first {it:#} non-omitted coefficients in {cmd:e(b)}.{p_end}

{phang}
{opt avc(#)} specifies the calibrated worst-case average correlation across observations (default {cmd:0.03}, strictly between {cmd:0.001} and {cmd:0.99}).{p_end}

{phang}
{opt cvs} stores and displays critical values for 32%, 10%, 5%, and 1% two-sided significance levels in {cmd:e(scpccvs)}.{p_end}

{phang}
{opt uncond} uses unconditional rather than regressor-conditioned critical values (required after {cmd:areg}).{p_end}

{phang}
{opt latlong} interprets {cmd:s_1} as latitude and {cmd:s_2} as longitude in degrees and uses spherical chordal distances.{p_end}

{phang}
{opt cluster(varname)} treats observations sharing a value of {it:varname} as a single spatial site (with constant site coordinates).{p_end}

{dlgtab:Geometry and Approximation}

{phang}
{opt method(string)}, {opt exact}, {opt approx} select exact eigendecomposition ({cmd:exact}), randomized Nyström landmark approximation ({cmd:nystrom}), or automatic selection ({cmd:auto}, switching to Nyström at 4,500+ sites).{p_end}

{phang}
{opt seed(#)}, {opt landmarks(#)}, {opt samples(#)}, {opt tolerance(#)} configure Nyström approximation and minimax critical-value optimization.{p_end}

{dlgtab:Output and Runtime}

{phang}
{opt weights} stores the site-level projection weight matrix {cmd:e(W)} and first representative row index per site {cmd:e(site_rows)}.{p_end}

{phang}
{opt backend(string)}, {opt nthreads(#)}, {opt capabilities} control or inspect the unified Rust plugin runtime.{p_end}


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. regress y x1 x2}{p_end}
{phang2}{cmd:. scpc}{p_end}
{phang2}{cmd:. scpc, latlong avc(0.05) cvs}{p_end}
{phang2}{cmd:. ivregress 2sls y (x1 = z1)}{p_end}
{phang2}{cmd:. scpc, latlong}{p_end}
{phang2}{cmd:. areg y x1, absorb(region)}{p_end}
{phang2}{cmd:. scpc, uncond}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:scpc} adds the following to {cmd:e()} while keeping {cmd:e(b)}, {cmd:e(V)}, and {cmd:e(sample)} intact (except {cmd:scpc, capabilities}, which only populates {cmd:r()}):

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(c0)}}calibrated kernel scale parameter at {opt avc}{p_end}
{synopt:{cmd:e(cmax)}}calibrated minimax scale parameter{p_end}
{synopt:{cmd:e(cv)}}reported 95% critical value (last printed coefficient){p_end}
{synopt:{cmd:e(q)}}number of SCPC low-frequency principal components selected{p_end}
{synopt:{cmd:e(tied)}}1 if boundary eigenvalues were tied, 0 otherwise{p_end}
{synopt:{cmd:e(n_obs)}}number of regression observations{p_end}
{synopt:{cmd:e(n_sites)}}number of spatial sites{p_end}
{synopt:{cmd:e(distance_scale)}}internal distance normalization scale{p_end}
{synopt:{cmd:e(avc)}}requested average correlation{p_end}
{synopt:{cmd:e(cv_tolerance)}}critical-value solver tolerance{p_end}
{synopt:{cmd:e(seed)}, {cmd:e(landmarks)}, {cmd:e(landmark_samples)}}Nyström settings (stored when {cmd:method(nystrom)} is used){p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(inference)}}{cmd:conditional_ols}, {cmd:conditional_2sls}, or {cmd:unconditional}{p_end}
{synopt:{cmd:e(method)}}{cmd:exact} or {cmd:nystrom}{p_end}
{synopt:{cmd:e(backend)}}{cmd:rust}{p_end}
{synopt:{cmd:e(specification)}}native specification version ({cmd:0.5.0}){p_end}
{synopt:{cmd:e(provider)}}linear algebra provider reported by the Rust core{p_end}
{synopt:{cmd:e(spatial_cluster)}}clustering variable name if clustered{p_end}
{synopt:{cmd:e(coordinate_units)}}coordinate distance units{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(scpcstats)}}table with columns {cmd:coef se t p ci_lower ci_upper}{p_end}
{synopt:{cmd:e(critical_95)}}95% critical values per reported coefficient{p_end}
{synopt:{cmd:e(scpccvs)}}critical values ({cmd:cv32 cv10 cv05 cv01}) when {opt cvs} is specified{p_end}
{synopt:{cmd:e(W)}, {cmd:e(site_rows)}}site projection weights and row indices when {opt weights} is specified{p_end}
{p2colreset}{...}


{marker authors}{...}
{title:Authors}

{pstd}
Paul David Boll, Sascha O. Becker, Ulrich K. Müller, Mark W. Watson, and Hans-Joachim Voth.{break}
Unified Rust architecture for Spatial SPUR.


{marker references}{...}
{title:References}

{marker mw2022}{...}
{phang}Müller, Ulrich K. and Mark W. Watson (2022). "Spatial Correlation Robust Inference." {it:Econometrica}, 90(6), 2901–2935.

{marker mw2023}{...}
{phang}Müller, Ulrich K. and Mark W. Watson (2023). "Spatial Correlation Robust Inference in Linear Regression and Panel Models." {it:Journal of Business & Economic Statistics}, 41(4), 1050–1064.
{p_end}
