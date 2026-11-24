{smcl}
{* *! version 0.6.0  24nov2026}{...}
{viewerjumpto "Syntax" "spurtest##syntax"}{...}
{viewerjumpto "Description" "spurtest##description"}{...}
{viewerjumpto "Options" "spurtest##options"}{...}
{viewerjumpto "Examples" "spurtest##examples"}{...}
{viewerjumpto "Stored results" "spurtest##results"}{...}
{viewerjumpto "References" "spurtest##references"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{bf:spurtest} {hline 2}}Individual Spatial Unit Root (SPUR) diagnostic tests{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:spurtest} {it:test} {varlist} {ifin} [{cmd:,} {it:options}]

{pstd}
where {it:test} is one of {cmd:i0}, {cmd:i1}, {cmd:i0resid}, or {cmd:i1resid}.
Levels diagnostics ({cmd:i0}, {cmd:i1}) take a single numeric variable; residual diagnostics ({cmd:i0resid}, {cmd:i1resid}) take {it:depvar} followed by one or more regressors.

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt q(#)}}number of low-frequency spatial directions ({cmd:2 <= q < n}); default is {cmd:15}{p_end}
{synopt :{opt nr:ep(#)}}number of Monte Carlo calibration draws; default is {cmd:100000}{p_end}
{synopt :{opt seed(#)}}random seed in {cmd:[0, 2147483647]}; default is {cmd:42}{p_end}
{synopt :{opt lat:long}}interpret {cmd:s_1}, {cmd:s_2} as latitude and longitude in degrees{p_end}
{synopt :{opt back:end(string)}}native backend ({cmd:rust} or {cmd:auto}){p_end}
{synopt :{opt nt:hreads(#)}}thread count override for the Rust backend{p_end}
{synopt :{opt drawfile(filename)}}precomputed little-endian float64 draws of dimension {cmd:2*q x nrep}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:spurtest} computes one of the four point-optimal spatial persistence tests from {help spurtest##mw2024:Müller and Watson (2024)} using spatial coordinates in {cmd:s_1}, {cmd:s_2}, ...:
{cmd:i0} tests the stationary {it:I}(0) null on a single series; {cmd:i1} tests the spatial unit-root {it:I}(1) null; {cmd:i0resid} and {cmd:i1resid} test the regression residual projection onto the exact Lévy–Brownian motion / stationary nuisance design.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. spurtest i0 y, q(15) latlong}{p_end}
{phang2}{cmd:. spurtest i1 y, q(15) latlong}{p_end}
{phang2}{cmd:. spurtest i0resid y x1 x2, q(15)}{p_end}
{phang2}{cmd:. spurtest i1resid y x1 x2, q(15)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spurtest} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(p)}}simulated p-value{p_end}
{synopt:{cmd:r(teststat)}}SPUR likelihood-ratio test statistic{p_end}
{synopt:{cmd:r(ha_param)}}selected normalized alternative scale parameter{p_end}
{synopt:{cmd:r(tied)}}1 if boundary eigenspaces were tied, 0 otherwise{p_end}
{synopt:{cmd:r(all_directions)}}1 if all positive directions were retained, 0 otherwise{p_end}
{synopt:{cmd:r(best_tuning_power)}}estimated power of the selected alternative grid point{p_end}
{synopt:{cmd:r(n_obs)}, {cmd:r(n_sites)}}number of observations and sites{p_end}
{synopt:{cmd:r(q)}, {cmd:r(nrep)}, {cmd:r(tuning_draws)}}test order and simulation draw counts{p_end}
{synopt:{cmd:r(seed)}}seed used when {opt drawfile()} is not specified{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(test)}}test name ({cmd:i0}, {cmd:i1}, {cmd:i0resid}, or {cmd:i1resid}){p_end}
{synopt:{cmd:r(normalization)}}{cmd:sample diameter}{p_end}
{synopt:{cmd:r(coordinate_units)}}{cmd:meters} or {cmd:coordinate units}{p_end}
{synopt:{cmd:r(backend)}, {cmd:r(specification)}, {cmd:r(provider)}}Rust runtime metadata{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:r(cv)}}3x1 matrix of critical values at the 1%, 5%, and 10% levels{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker mw2024}{...}
{phang}Müller, Ulrich K. and Mark W. Watson (2024). "Spatial Unit Roots and Spurious Regression." {it:Econometrica}, 92(5), 1661–1695.
{p_end}
