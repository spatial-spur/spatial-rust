clear all
set more off
args package plus
sysdir set PLUS "`plus'"
net install spatial_spur, from("`package'") replace
import delimited using "common.csv", asdouble clear
gen double s_1 = lat
gen double s_2 = lon
spur y x x2, q(6) nrep(400) drawfile("diagnostics.bin") latlong cvs backend(rust) nthreads(1)
assert "`r(backend)'" == "rust"

tempname stats levels adjusted cvlevels cvadjusted cvs_levels cvs_adjusted diagnostics cuts handles
matrix `levels' = r(levels_scpcstats)
matrix `adjusted' = r(transformed_scpcstats)
matrix `cvlevels' = r(levels_critical_95)
matrix `cvadjusted' = r(transformed_critical_95)
matrix `cvs_levels' = r(levels_scpccvs)
matrix `cvs_adjusted' = r(transformed_scpccvs)
matrix `diagnostics' = r(diagnostics)
matrix `cuts' = r(diagnostic_cutoffs)
local q_levels = r(levels_q)
local q_adjusted = r(transformed_q)
local n_common = r(n_obs)
postfile `stats' str20 case str9 stage str10 term int n int q double(coef se t p lo hi cv) using scpc-stata.dta
postfile `handles' str9 stage str10 term double(cv32 cv10 cv05 cv01) using cvs-stata.dta
forvalues j = 1/3 {
    local term : word `j' of x x2 _cons
    post `stats' ("common") ("levels") ("`term'") (`n_common') (`q_levels') ///
       (`levels'[`j',1]) (`levels'[`j',2]) (`levels'[`j',3]) (`levels'[`j',4]) ///
       (`levels'[`j',5]) (`levels'[`j',6]) (`cvlevels'[`j',1])
    post `stats' ("common") ("adjusted") ("`term'") (`n_common') (`q_adjusted') ///
       (`adjusted'[`j',1]) (`adjusted'[`j',2]) (`adjusted'[`j',3]) (`adjusted'[`j',4]) ///
       (`adjusted'[`j',5]) (`adjusted'[`j',6]) (`cvadjusted'[`j',1])
    post `handles' ("levels") ("`term'") ///
       (`cvs_levels'[`j',1]) (`cvs_levels'[`j',2]) (`cvs_levels'[`j',3]) (`cvs_levels'[`j',4])
    post `handles' ("adjusted") ("`term'") ///
       (`cvs_adjusted'[`j',1]) (`cvs_adjusted'[`j',2]) (`cvs_adjusted'[`j',3]) (`cvs_adjusted'[`j',4])
}
postclose `handles'
postfile `handles' str8 test double(statistic p ha cv01 cv05 cv10) using diagnostics-stata.dta
local j = 0
foreach test in i0 i1 i0resid i1resid {
    local ++j
    post `handles' ("`test'") (`diagnostics'[`j',1]) (`diagnostics'[`j',2]) (`diagnostics'[`j',3]) ///
      (`cuts'[`j',1]) (`cuts'[`j',2]) (`cuts'[`j',3])
}
postclose `handles'

spurhalflife y, q(6) nrep(400) drawfile("halflife.bin") latlong normdist backend(rust) nthreads(1)
local components = r(n_components)
postfile `handles' double(lower upper) using halflife-stata.dta
if `components' > 0 {
    matrix bounds = r(components)
    forvalues i = 1/`components' {
        post `handles' (bounds[`i',1]) (bounds[`i',2])
    }
}
postclose `handles'

spurhalflife y, q(6) nrep(400) drawfile("halflife.bin") latlong backend(rust) nthreads(1)
local components_raw = r(n_components)
postfile `handles' double(lower upper) using halflife_raw-stata.dta
if `components_raw' > 0 {
    matrix bounds = r(components)
    forvalues i = 1/`components_raw' {
        post `handles' (bounds[`i',1]) (bounds[`i',2])
    }
}
postclose `handles'

preserve
spurtransform y, prefix(h_) transformation(nn) latlong backend(rust) nthreads(1) replace
keep h_y
save tr_nn_stata.dta, replace
restore
preserve
spurtransform y, prefix(h_) transformation(iso) radius(500000) latlong backend(rust) nthreads(1) replace
keep h_y
save tr_iso_stata.dta, replace
restore

import delimited using "clustered.csv", asdouble clear
gen double s_1 = lat
gen double s_2 = lon
quietly regress y x
quietly scpc, cluster(cl) latlong backend(rust) nthreads(1)
matrix `levels' = e(scpcstats)
matrix `cvlevels' = e(critical_95)
forvalues j = 1/2 {
    local term : word `j' of x _cons
    post `stats' ("clustered") ("levels") ("`term'") (e(n_obs)) (e(q)) ///
      (`levels'[`j',1]) (`levels'[`j',2]) (`levels'[`j',3]) (`levels'[`j',4]) ///
      (`levels'[`j',5]) (`levels'[`j',6]) (`cvlevels'[`j',1])
}
spurtransform y, prefix(h_) transformation(cluster) clustvar(cl) backend(rust) nthreads(1)
keep h_y
save tr_cl_stata.dta, replace
import delimited using "common.csv", asdouble clear
gen double s_1 = lat
gen double s_2 = lon

capture confirm file "chetty-cases.csv"
if !_rc {
    import delimited using "chetty-cases.csv", asdouble clear
    gen double s_1 = lat
    gen double s_2 = lon
    levelsof term, local(keys) clean
    foreach key of local keys {
        quietly regress y x if term == "`key'", nocons
        quietly scpc, latlong nthreads(1)
        matrix `levels' = e(scpcstats)
        matrix `cvlevels' = e(critical_95)
        post `stats' ("`key'") ("levels") ("x") (e(n_obs)) (e(q)) ///
          (`levels'[1,1]) (`levels'[1,2]) (`levels'[1,3]) (`levels'[1,4]) ///
          (`levels'[1,5]) (`levels'[1,6]) (`cvlevels'[1,1])
        quietly spurtransform y x if term == "`key'", prefix(h_) replace latlong nthreads(1)
        quietly regress h_y h_x if term == "`key'", nocons
        quietly scpc, latlong nthreads(1)
        matrix `adjusted' = e(scpcstats)
        matrix `cvadjusted' = e(critical_95)
        post `stats' ("`key'") ("adjusted") ("x") (e(n_obs)) (e(q)) ///
          (`adjusted'[1,1]) (`adjusted'[1,2]) (`adjusted'[1,3]) (`adjusted'[1,4]) ///
          (`adjusted'[1,5]) (`adjusted'[1,6]) (`cvadjusted'[1,1])
    }
}
postclose `stats'
di "Stata parity completed"
