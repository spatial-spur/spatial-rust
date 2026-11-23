{smcl}
{title:Spatial unit root diagnostics}

{p 4 8 2}{cmd:spurtest} {it:i0|i1|i0resid|i1resid y [xvars]}{cmd:,} [{cmd:latlong} {cmd:q(}{it:#}{cmd:)} {cmd:nrep(}{it:#}{cmd:)}]{p_end}

{p 4 4 2}Use {cmd:i0} or {cmd:i1} for levels or their residual variants for an OLS design with an intercept.{p_end}

{p 4 4 2}{cmd:drawfile()} accepts a row-major little-endian double matrix with {cmd:2*q} rows and {cmd:nrep} columns, inference rows first and independent tuning rows second.{p_end}
