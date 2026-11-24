{smcl}
{* *! version 0.6.0  24nov2026}{...}
{viewerjumpto "Syntax" "spur##syntax"}{...}
{viewerjumpto "Description" "spur##description"}{...}
{viewerjumpto "Options" "spur##options"}{...}
{viewerjumpto "Examples" "spur##examples"}{...}
{viewerjumpto "Stored results" "spur##results"}{...}
{viewerjumpto "References" "spur##references"}{...}
{title:Title}

{p2colset 5 14 16 2}{...}
{p2col :{bf:spur} {hline 2}}Combined Spatial Unit Root (SPUR) diagnostics, LBM-GLS transformation, and SCPC inference{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:spur} {it:depvar indepvars} {ifin} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:SPUR Diagnostics & Transformation}
{synopt :{opt q(#)}}number of low-frequency spatial directions ({cmd:>= 2}); default is {cmd:15}{p_end}
{synopt :{opt nr:ep(#)}}Monte Carlo draws for critical values and tuning; default is {cmd:100000}{p_end}
{synopt :{opt seed(#)}}random seed in {cmd:[0, 2147483647]}; default is {cmd:42}{p_end}
{synopt :{opt prefix(name)}}prefix for generated LBM-GLS transformed variables; default is {cmd:h_}{p_end}
{synopt :{opt replace}}overwrite existing variables sharing {opt prefix()}{p_end}
{synopt :{opt drawfile(filename)}}precomputed little-endian float64 draws of dimension {cmd:2*q x nrep}{p_end}

{syntab:SCPC & Runtime}
{synopt :{opt lat:long}}interpret {cmd:s_1}, {cmd:s_2} as latitude and longitude in degrees{p_end}
{synopt :{opt avc(#)}}SCPC average pairwise correlation in {cmd:(0.001, 0.99)}; default is {cmd:0.03}{p_end}
{synopt :{opt unc:ond}}use unconditional SCPC critical values{p_end}
{synopt :{opt cvs}}store full SCPC critical value tables ({cmd:32%, 10%, 5%, 1%}){p_end}
{synopt :{opt method(string)}}SCPC geometry method: {cmd:auto}, {cmd:exact}, or {cmd:nystrom}{p_end}
{synopt :{opt back:end(string)}}native backend ({cmd:rust} or {cmd:auto}){p_end}
{synopt :{opt nt:hreads(#)}}thread count override for the Rust backend{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:spur} implements the full spatial unit root diagnostic and correction workflow of {help spur##mw2024:Müller and Watson (2024)} and {help spur##bbv2026:Becker, Boll, and Voth (2026)}.
Given coordinates in {cmd:s_1}, {cmd:s_2}, ..., {cmd:spur}: (1) runs the four SPUR diagnostics ({cmd:i0}, {cmd:i1}, {cmd:i0resid}, {cmd:i1resid}), (2) creates LBM-GLS spatially demeaned/transformed variables prefixed with {opt prefix()}, and (3) estimates SCPC inference for both the levels regression and the LBM-GLS transformed regression.

{pstd}
Note that LBM-GLS transformed variables are weighted-centered rather than untransformed, and {opt q()} must satisfy {cmd:2 <= q <= n - p - 1}.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. spur y x1 x2}{p_end}
{phang2}{cmd:. spur y x1 x2, latlong q(15) avc(0.03) prefix(h_) replace cvs}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spur} stores the following in {cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:r(i0_p)}, {cmd:r(i0_teststat)}, {cmd:r(i0_ha_param)}}levels I(0) p-value, statistic, and alternative tuning parameter{p_end}
{synopt:{cmd:r(i1_p)}, {cmd:r(i1_teststat)}, {cmd:r(i1_ha_param)}}levels I(1) p-value, statistic, and alternative tuning parameter{p_end}
{synopt:{cmd:r(i0resid_p)}, {cmd:r(i0resid_teststat)}, {cmd:r(i0resid_ha_param)}}residual I(0) p-value, statistic, and alternative parameter{p_end}
{synopt:{cmd:r(i1resid_p)}, {cmd:r(i1resid_teststat)}, {cmd:r(i1resid_ha_param)}}residual I(1) p-value, statistic, and alternative parameter{p_end}
{synopt:{cmd:r(q)}, {cmd:r(nrep)}, {cmd:r(tuning_draws)}}SPUR diagnostic order and Monte Carlo draw counts{p_end}
{synopt:{cmd:r(n_obs)}, {cmd:r(n_sites)}}observation and spatial site counts{p_end}
{synopt:{cmd:r(levels_q)}, {cmd:r(transformed_q)}}SCPC {it:q} for levels and transformed regressions{p_end}
{synopt:{cmd:r(seed)}}Monte Carlo seed when {opt drawfile()} is not given{p_end}

{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:r(outputs)}}names of generated LBM-GLS variables{p_end}
{synopt:{cmd:r(backend)}, {cmd:r(specification)}, {cmd:r(provider)}}Rust backend metadata{p_end}

{p2col 5 28 32 2: Matrices}{p_end}
{synopt:{cmd:r(diagnostics)}}4x3 matrix ({cmd:statistic p ha_param}) for {cmd:i0 i1 i0resid i1resid}{p_end}
{synopt:{cmd:r(diagnostic_info)}}4x3 matrix ({cmd:tied all_directions best_tuning_power}){p_end}
{synopt:{cmd:r(diagnostic_cutoffs)}}4x3 matrix of 1%, 5%, 10% critical values{p_end}
{synopt:{cmd:r(levels_scpcstats)}, {cmd:r(transformed_scpcstats)}}SCPC regression tables before and after LBM-GLS{p_end}
{synopt:{cmd:r(levels_critical_95)}, {cmd:r(transformed_critical_95)}}95% SCPC critical values{p_end}
{synopt:{cmd:r(levels_scpccvs)}, {cmd:r(transformed_scpccvs)}}additional critical values when {opt cvs} is specified{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker bbv2026}{...}
{phang}Becker, Sascha O., P. David Boll, and Hans-Joachim Voth (2026). "Testing and Correcting for Spatial Unit Roots in Regression Analysis." {it:The Stata Journal}, 26(2), 177–202.

{marker mw2024}{...}
{phang}Müller, Ulrich K. and Mark W. Watson (2024). "Spatial Unit Roots and Spurious Regression." {it:Econometrica}, 92(5), 1661–1695.
{p_end}
