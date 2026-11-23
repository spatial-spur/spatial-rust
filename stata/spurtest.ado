*! Spatial unit root diagnostic
program spurtest, rclass sortpreserve
    version 16.0
    gettoken test 0 : 0, parse(" ,")
    syntax varlist(numeric) [if] [in], [Q(integer 15) NRep(integer 100000) SEED(real 42) ///
        LATlong BACKend(string) NThreads(integer -1) DRAWfile(string)]
    if !inlist("`test'", "i0", "i1", "i0resid", "i1resid") {
        di as error "Choose i0, i1, i0resid, or i1resid"
        exit 198
    }
    local residual = inlist("`test'", "i0resid", "i1resid")
    gettoken y x : varlist
    local p : word count `x'
    if !`residual' & `p' {
        di as error "A levels diagnostic takes one response variable"
        exit 198
    }
    marksample touse
    _spatial_spur load, backend(`"`backend'"') nthreads(`nthreads')
    local provider "`r(provider)'"
    local specification "`r(specification)'"
    _spatial_spur coords, sample(`touse') `latlong'
    local coords `r(coords)'
    local d = r(d)
    local n = r(n)
    if `q' < 2 | `q' >= `n' | `nrep' < 1 {
        di as error "q() must be at least two and smaller than the sample; nrep() must be positive"
        exit 198
    }
    local use_file `"`drawfile'"'
    if `"`use_file'"' == "" {
        tempfile normals
        _spatial_spur draws, q(`=2*`q'') nrep(`nrep') seed(`seed') file(`"`normals'"')
        local use_file `"`normals'"'
    }
    local native = "spur_" + subinstr("`test'", "resid", "_resid", .)
    local ll = "`latlong'" != ""
    tempname critical statistic pvalue alternative info
    matrix `critical' = J(3,1,.)
    matrix `info' = J(1,3,.)
    _spatial_spur invoke, vars("`coords' `y' `x'") sample(`touse') ///
        command(`"`native' `d' `p' `residual' `ll' `q' `nrep' "`nthreads'" "`use_file'" `critical' `statistic' `pvalue' `alternative' `info'"')
    di as text "SPUR `test' p-value " as result %8.5f scalar(`pvalue')
    return scalar p = scalar(`pvalue')
    return scalar teststat = scalar(`statistic')
    return scalar ha_param = scalar(`alternative')
    return matrix cv = `critical'
    return scalar tied = `info'[1,1]
    return scalar all_directions = `info'[1,2]
    return scalar best_tuning_power = `info'[1,3]
    return scalar n_obs = `n'
    return scalar n_sites = `n'
    return scalar q = `q'
    return scalar nrep = `nrep'
    return scalar tuning_draws = `nrep'
    if `"`drawfile'"' == "" return scalar seed = `seed'
    return local test "`test'"
    return local normalization "sample diameter"
    return local coordinate_units = cond(`ll', "meters", "coordinate units")
    return local backend "rust"
    return local specification "`specification'"
    return local provider "`provider'"
end
