{smcl}
{title:Inference for fitted spatial regressions}

{p 4 8 2}{cmd:scpc}{cmd:,} [{cmd:latlong} {cmd:uncond} {cmd:cvs} {cmd:weights} {cmd:method(}{it:choice}{cmd:)} {cmd:cluster(}{it:varname}{cmd:)}]{p_end}

{p 4 4 2}Run after unweighted {cmd:regress}, ordinary {cmd:ivregress 2sls} without absorption, or one-way {cmd:areg} with {cmd:uncond}. The {cmd:areg} results report slopes.{p_end}

{p 4 4 2}{cmd:s_1}, {cmd:s_2}, and successive numbered variables hold the spatial coordinates. For {cmd:latlong}, latitude is first.{p_end}

{p 4 4 2}{cmd:method()} accepts {cmd:auto}, {cmd:exact}, or {cmd:nystrom}. Nyström accepts {cmd:seed()}, {cmd:landmarks()}, and {cmd:samples()}.{p_end}

{p 4 4 2}{cmd:e(critical_95)} contains the effective critical value for each reported coefficient. {cmd:e(cv)} is the location-only cutoff. {cmd:weights} saves {cmd:e(W)} and each row's representative observation number in {cmd:e(site_rows)}.{p_end}
