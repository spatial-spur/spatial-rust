clear all
set more off
args package plus
sysdir set PLUS "`plus'"
net install spatial_spur, from("`package'") replace

mata:
void spatial_spur_test_weights() {
  real matrix X, W, V
  real colvector y
  real rowvector s
  X = st_data(., "x"), J(st_nobs(),1,1)
  y = st_data(., st_local("epsilon"))
  W = st_matrix("W")
  W = W[,2..cols(W)]
  V = W' * ((X :* y) * invsym(quadcross(X,X)))
  s = sqrt(st_nobs()/st_numscalar("n_q")) :* sqrt(colsum(V:^2))
  assert(max(abs(s' - st_matrix("table")[,2])) < 1e-10)
}
end
set obs 96
set seed 9082
gen double s_1 = runiform()
gen double s_2 = runiform()
gen int g = mod(_n, 8)
gen int h = ceil(_n/16)
gen double z = rnormal()
gen double x = z + rnormal()
gen double y = x + rnormal() + g + h
gen byte h_yy = 7
set varabbrev on

local state `"`c(rngstate)'"'
spurtest i0 y, q(4) nrep(200) backend(rust)
assert "`r(backend)'" == "rust"
assert `"`c(rngstate)'"' == `"`state'"'
gen double const = 1
capture quietly spurtest i0 const, q(4) nrep(200)
assert _rc == 198
assert `"`c(rngstate)'"' == `"`state'"'
drop const
spurhalflife y, q(4) nrep(200) normdist backend(rust) nthreads(1)
assert r(normdist) == 1
assert r(n_components) > 0
assert "`r(distance_units)'" == "coordinate units"
assert `"`c(rngstate)'"' == `"`state'"'

spurtransform y, prefix(h_) replace transformation(nn) nthreads(1)
assert h_yy == 7
assert !missing(h_y)
quietly datasignature
local before = r(datasignature)
capture quietly spur y x, q(4) nrep(200) replace avc(.999) nthreads(1)
assert _rc == 198
quietly datasignature
assert "`r(datasignature)'" == "`before'"
assert `"`c(rngstate)'"' == `"`state'"'
gen double abcdefghijklmnopqrstuvwxyzabcdef = y
quietly datasignature
local before = r(datasignature)
capture quietly spurtransform y abcdefghijklmnopqrstuvwxyzabcdef, prefix(h_) replace transformation(nn)
assert _rc == 198
quietly datasignature
assert "`r(datasignature)'" == "`before'"

quietly regress y x
spur y x, q(4) nrep(200) replace cvs nthreads(1)
assert "`e(cmd)'" == "regress"
assert "`r(backend)'" == "rust"
assert colsof(r(levels_scpcstats)) == 6
assert rowsof(r(levels_scpcstats)) == 2
assert rowsof(r(transformed_scpcstats)) == 1
assert !missing(h_x)
assert h_yy == 7
assert `"`c(rngstate)'"' == `"`state'"'

quietly regress y x
preserve
quietly replace y = y + 1
capture quietly scpc
assert _rc == 198
restore

quietly regress y x
preserve
tempvar u_orig
quietly predict double `u_orig', residuals
quietly replace y = _b[_cons] + _b[x] * x + `u_orig'[cond(mod(_n,2)==1, _n+1, _n-1)]
capture quietly scpc
assert _rc != 0
restore

quietly areg y x, absorb(g)
capture quietly scpc
assert _rc == 198
quietly scpc, uncond nthreads(1)
matrix slopes = e(scpcstats)
quietly regress y x i.g
quietly scpc, uncond k(1) nthreads(1)
matrix dummy = e(scpcstats)
assert mreldif(dummy, slopes) < 1e-11

capture quietly areg y x, absorb(g h)
if !_rc {
    capture quietly scpc, uncond
    assert _rc == 198
}
capture quietly ivregress 2sls y (x=z), absorb(g)
if !_rc {
    capture quietly scpc, uncond
    assert _rc == 198
}
quietly ivregress 2sls y (x=z)
quietly scpc, cvs nthreads(1)
assert "`e(inference)'" == "conditional_2sls"
matrix iv = e(scpcstats)
assert abs(iv[1,1] - _b[x]) < 1e-14

quietly regress y x
quietly scpc, method(nystrom) landmarks(52) samples(4) seed(7) weights nthreads(1)
assert "`e(method)'" == "nystrom"
assert e(seed) == 7
assert rowsof(e(W)) == _N
assert colsof(e(W)) == e(q) + 1
matrix W = e(W)
matrix b = e(b)
matrix table = e(scpcstats)
scalar n_q = e(q)
tempvar epsilon
quietly predict double `epsilon', residuals
mata: spatial_spur_test_weights()
drop `epsilon'

preserve
clear
set obs 17
gen double s_1 = 0
gen double s_2 = -120 + (_n-1)/4
gen double y = _n
spurtransform y, prefix(h_) transformation(nn) latlong nthreads(1)
assert h_y == 0 if inrange(_n,2,16)
assert h_y == -1 if _n == 1
assert h_y == 1 if _n == 17
replace y = 1e12
spurtransform y, prefix(h_) transformation(iso) latlong radius(1000000) nthreads(1) replace
assert h_y == 0
restore

scpc, capabilities
assert "`r(specification)'" == "0.5.0"
assert "`r(backend)'" == "rust"

spurhalflife y, q(4) nrep(200) normdist backend(rust) nthreads(1)
matrix comp_norm = r(components)
local max_d = r(max_dist)
spurhalflife y, q(4) nrep(200) backend(rust) nthreads(1)
assert r(normdist) == 0
assert r(n_components) > 0
matrix comp_raw = r(components)
matrix comp_scaled = comp_norm * `max_d'
assert mreldif(comp_raw, comp_scaled) < 1e-12

mata: mata clear
mata: _sp_fd = 12345
spurtest i0 y, q(4) nrep(200) backend(rust) nthreads(1)
spurhalflife y, q(4) nrep(200) backend(rust) nthreads(1)
quietly regress y x
scpc, cvs method(nystrom) landmarks(52) samples(4) seed(7) weights nthreads(1)
mata: assert(_sp_fd == 12345)
mata: mata drop _sp_fd
assert colsof(e(scpccvs)) == 4
assert e(seed) == 7
scpc, method(exact) nthreads(1)
assert "`e(cmd)'" == "regress"
capture confirm matrix e(scpccvs)
assert _rc != 0
capture confirm matrix e(W)
assert _rc != 0
capture confirm scalar e(seed)
assert _rc != 0

preserve
gen str4 sg = string(g)
quietly areg y x, absorb(g)
quietly scpc, uncond nthreads(1)
matrix slopes_num = e(scpcstats)
quietly areg y x, absorb(sg)
quietly scpc, uncond nthreads(1)
matrix slopes_str = e(scpcstats)
assert mreldif(slopes_num, slopes_str) < 1e-12
restore

preserve
replace g = . in 1/3
capture quietly spurtransform y in 1/3, prefix(cl_) transformation(cluster) clustvar(g) replace
assert _rc == 2001
restore

capture quietly spur y x, q(4) nrep(200) replace method(bogus)
assert _rc == 198

di "Stata installed integration checks passed"
