{smcl}
{title:Run the combined SPUR workflow}

{p 4 8 2}{cmd:spur} {it:y xvars} [{it:if}] [{it:in}]{cmd:,} [{cmd:latlong} {cmd:q(}{it:#}{cmd:)} {cmd:nrep(}{it:#}{cmd:)} {cmd:replace}]{p_end}

{p 4 4 2}Runs four diagnostics, applies LBM-GLS, and returns levels and transformed SCPC tables. Outputs with prefix {cmd:h_} are committed after the two regressions succeed.{p_end}

{p 4 4 2}The seed leaves the caller's RNG state unchanged. See {help spurtest} for the common-draw format.{p_end}
