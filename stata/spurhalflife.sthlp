{smcl}
{* *! version 0.6.0  24nov2026}{...}
{viewerjumpto "Syntax" "spurhalflife##syntax"}{...}
{viewerjumpto "Description" "spurhalflife##description"}{...}
{viewerjumpto "Options" "spurhalflife##options"}{...}
{viewerjumpto "Examples" "spurhalflife##examples"}{...}
{viewerjumpto "Stored results" "spurhalflife##results"}{...}
{viewerjumpto "References" "spurhalflife##references"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{bf:spurhalflife} {hline 2}}Spatial half-life confidence set by SPUR test inversion{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:spurhalflife} {varname} {ifin} [{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt q(#)}}number of low-frequency spatial directions ({cmd:2 <= q <= 100}); default is {cmd:15}{p_end}
{synopt :{opt nr:ep(#)}}number of Monte Carlo draws; default is {cmd:100000}{p_end}
{synopt :{opt seed(#)}}random seed in {cmd:[0, 2147483647]}; default is {cmd:42}{p_end}
{synopt :{opt level(#)}}confidence level in {cmd:(0, 100)}; default is {cmd:95}{p_end}
{synopt :{opt lat:long}}interpret {cmd:s_1}, {cmd:s_2} as latitude and longitude in degrees{p_end}
{synopt :{opt norm:dist}}report endpoints divided by sample diameter instead of raw coordinate units / meters{p_end}
{synopt :{opt back:end(string)}}native backend ({cmd:rust} or {cmd:auto}){p_end}
{synopt :{opt nt:hreads(#)}}thread count override for the Rust backend{p_end}
{synopt :{opt drawfile(filename)}}precomputed little-endian float64 draws of dimension {cmd:q x nrep}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:spurhalflife} computes the continuous-domain SPUR confidence set for the exponential spatial decay half-life (including boundary limits at 0 and infinity and reporting all accepted connected components in {cmd:r(components)}) using coordinates in {cmd:s_1}, {cmd:s_2}, ....
By default, endpoints are expressed in coordinate units (or meters with {opt latlong}); specify {opt normdist} to express them as multiples of the sample diameter {cmd:r(max_dist)}.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. spurhalflife y}{p_end}
{phang2}{cmd:. spurhalflife y, q(15) level(95) latlong normdist}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spurhalflife} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:r(ci_l)}, {cmd:r(ci_u)}}outer convex-hull lower and upper endpoints ({cmd:.} for infinity){p_end}
{synopt:{cmd:r(max_dist)}}maximum pairwise sample distance (diameter){p_end}
{synopt:{cmd:r(n_components)}}number of connected intervals in the confidence set{p_end}
{synopt:{cmd:r(empty)}, {cmd:r(unbounded)}}indicators for an empty or right-unbounded confidence set{p_end}
{synopt:{cmd:r(tied)}, {cmd:r(all_directions)}}eigenspace tie and direction retention indicators{p_end}
{synopt:{cmd:r(n_obs)}, {cmd:r(n_sites)}}observation and site counts{p_end}
{synopt:{cmd:r(q)}, {cmd:r(nrep)}, {cmd:r(level)}, {cmd:r(normdist)}}specification parameters{p_end}
{synopt:{cmd:r(seed)}}seed used when {opt drawfile()} is not given{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:r(distance_units)}}{cmd:meters} or {cmd:coordinate units}{p_end}
{synopt:{cmd:r(normalization)}}{cmd:sample diameter}{p_end}
{synopt:{cmd:r(backend)}, {cmd:r(specification)}, {cmd:r(provider)}}Rust runtime metadata{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:r(components)}}{it:n_components} x 2 matrix of component bounds ({cmd:lower upper}) when non-empty{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker mw2024}{...}
{phang}Müller, Ulrich K. and Mark W. Watson (2024). "Spatial Unit Roots and Spurious Regression." {it:Econometrica}, 92(5), 1661–1695.
{p_end}
