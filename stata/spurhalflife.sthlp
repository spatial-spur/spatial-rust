{smcl}
{title:Spatial half-life confidence sets}

{p 4 8 2}{cmd:spurhalflife} {it:varname}{cmd:,} [{cmd:latlong} {cmd:normdist} {cmd:q(}{it:#}{cmd:)} {cmd:nrep(}{it:#}{cmd:)}]{p_end}

{p 4 4 2}{cmd:r(components)} contains detected confidence components. A missing upper endpoint in a component represents infinity.{p_end}

{p 4 4 2}Geographic output uses meters. {cmd:normdist} uses sample diameters. A {cmd:drawfile()} contains {cmd:q} rows of little-endian Gaussian doubles in row-major order.{p_end}
