{smcl}
{* *! version 0.6.0  24nov2026}{...}
{viewerjumpto "Syntax" "spurtransform##syntax"}{...}
{viewerjumpto "Description" "spurtransform##description"}{...}
{viewerjumpto "Options" "spurtransform##options"}{...}
{viewerjumpto "Examples" "spurtransform##examples"}{...}
{viewerjumpto "Stored results" "spurtransform##results"}{...}
{viewerjumpto "References" "spurtransform##references"}{...}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{bf:spurtransform} {hline 2}}Spatial nuisance and unit-root transformations (LBM-GLS, NN, ISO, cluster){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 20 2}
{cmd:spurtransform} {varlist} {ifin} {cmd:,} {opt pre:fix(name)} [{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt pre:fix(name)}}required prefix for generated transformed variables{p_end}
{synopt :{opt trans:formation(string)}}transformation type: {cmd:lbmgls} (default), {cmd:nn}, {cmd:iso}, or {cmd:cluster}{p_end}
{synopt :{opt radius(#)}}positive neighborhood radius for {cmd:transformation(iso)}{p_end}
{synopt :{opt clustvar(varname)}}clustering variable for {cmd:transformation(cluster)}{p_end}
{synopt :{opt lat:long}}interpret {cmd:s_1}, {cmd:s_2} as latitude and longitude in degrees (meters for {opt radius()}){p_end}
{synopt :{opt replace}}overwrite existing variables sharing {opt prefix()}{p_end}
{synopt :{opt sep:arately}}evaluate complete observations separately for each variable in {varlist}{p_end}
{synopt :{opt back:end(string)}}native backend ({cmd:rust} or {cmd:auto}){p_end}
{synopt :{opt nt:hreads(#)}}thread count override for the Rust backend{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:spurtransform} computes spatially transformed versions of variables in {varlist}, reading coordinates from {cmd:s_1}, {cmd:s_2}, ... (except for {cmd:transformation(cluster)}) and staging outputs atomically with {opt prefix()}:

{phang2}{cmd:lbmgls}: Lévy–Brownian motion generalized least squares spatial whitening/demeaning.{p_end}
{phang2}{cmd:nn}: Nearest-neighbor demeaning (for tie-free neighbor groups).{p_end}
{phang2}{cmd:iso}: Uniform-radius spatial neighborhood demeaning using {opt radius()}.{p_end}
{phang2}{cmd:cluster}: Within-cluster demeaning using {opt clustvar()}.{p_end}


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. spurtransform y x1 x2, prefix(h_) replace}{p_end}
{phang2}{cmd:. spurtransform y, prefix(nn_) transformation(nn) latlong}{p_end}
{phang2}{cmd:. spurtransform y, prefix(iso_) transformation(iso) latlong radius(50000)}{p_end}
{phang2}{cmd:. spurtransform y, prefix(cl_) transformation(cluster) clustvar(region)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spurtransform} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:r(outputs)}}list of created output variable names{p_end}
{synopt:{cmd:r(transformation)}}transformation applied ({cmd:lbmgls}, {cmd:nn}, {cmd:iso}, or {cmd:cluster}){p_end}
{synopt:{cmd:r(backend)}, {cmd:r(specification)}, {cmd:r(provider)}}Rust runtime metadata{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker bbv2026}{...}
{phang}Becker, Sascha O., P. David Boll, and Hans-Joachim Voth (2026). "Testing and Correcting for Spatial Unit Roots in Regression Analysis." {it:The Stata Journal}, 26(2), 177–202.

{marker mw2024}{...}
{phang}Müller, Ulrich K. and Mark W. Watson (2024). "Spatial Unit Roots and Spurious Regression." {it:Econometrica}, 92(5), 1661–1695.
{p_end}
