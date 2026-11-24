*! Combined SPUR diagnostics and SCPC
program spur, rclass sortpreserve
    version 16.0
    syntax varlist(numeric min=2) [if] [in], [Q(integer 15) NRep(integer 100000) ///
        SEED(real 42) LATlong AVC(real .03) UNCond CVS REPLACE BACKend(string) ///
        NThreads(integer -1) DRAWfile(string) METHOD(string) PREFIX(name)]
    if "`prefix'" == "" local prefix h_
    if `q' < 2 | `nrep' < 1 {
        di as error "q() must be at least two and nrep() must be positive"
        exit 198
    }
    if !(`avc' > 0.001 & `avc' < 0.99) {
        di as error "avc() must be strictly between 0.001 and 0.99"
        exit 198
    }
    local method = lower(trim(`"`method'"'))
    if !inlist("`method'", "", "auto", "exact", "nystrom", "approx") {
        di as error "Use method(auto), method(exact), or method(nystrom)"
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

    gettoken y x : varlist
    local p : word count `x'
    if `q' > `n' - `p' - 1 {
        di as error "q() exceeds the residual degrees of freedom"
        exit 198
    }
    local outputs
    local source
    foreach v of local varlist {
        local outputs `outputs' `prefix'`v'
        tempvar result
        quietly gen double `result' = .
        local source `source' `result'
    }
    _spatial_spur validate, outputs("`outputs'") `replace'
    local use_file `"`drawfile'"'
    if `"`use_file'"' == "" {
        tempfile normals
        _spatial_spur draws, q(`=2*`q'') nrep(`nrep') seed(`seed') file(`"`normals'"')
        local use_file `"`normals'"'
    }
    local ll = "`latlong'" != ""
    tempname diagnostics cuts info prior
    matrix `diagnostics' = J(4,3,.)
    matrix `cuts' = J(4,3,.)
    matrix `info' = J(4,3,.)
    _spatial_spur invoke, vars("`coords' `varlist' `source'") sample(`touse') ///
        command(`"spur_pipeline `d' `p' 1 `ll' `q' `nrep' "`nthreads'" "`use_file'" `diagnostics' `cuts' `info'"')

    _estimates hold `prior', restore nullok
    local scpc_opts avc(`avc') `latlong' `uncond' `cvs' backend(`"`backend'"') nthreads(`nthreads')
    if "`method'" != "" local scpc_opts `scpc_opts' method(`method')
    quietly regress `varlist' if `touse'
    quietly scpc, `scpc_opts'
    tempname levels adjusted cv_levels cv_adjusted extra_levels extra_adjusted
    matrix `levels' = e(scpcstats)
    matrix `cv_levels' = e(critical_95)
    if "`cvs'" != "" matrix `extra_levels' = e(scpccvs)
    local levels_q = e(q)
    quietly regress `source' if `touse'
    quietly scpc, `scpc_opts'
    matrix `adjusted' = e(scpcstats)
    matrix `cv_adjusted' = e(critical_95)
    if "`cvs'" != "" matrix `extra_adjusted' = e(scpccvs)
    local transformed_q = e(q)
    local levels_names : rownames `levels'
    if rowsof(`adjusted') != `: word count `levels_names'' {
        _estimates unhold `prior'
        di as error "Transformed regressors became collinear after spatial transformation"
        exit 198
    }
    matrix rownames `adjusted' = `levels_names'
    matrix rownames `cv_adjusted' = `levels_names'
    if "`cvs'" != "" matrix rownames `extra_adjusted' = `levels_names'

    matrix rownames `diagnostics' = i0 i1 i0resid i1resid
    matrix colnames `diagnostics' = statistic p ha_param
    matrix rownames `info' = i0 i1 i0resid i1resid
    matrix colnames `info' = tied all_directions best_tuning_power
    matrix rownames `cuts' = i0 i1 i0resid i1resid
    matrix colnames `cuts' = cv01 cv05 cv10
    _estimates unhold `prior'
    _spatial_spur commit, outputs("`outputs'") source("`source'") `replace'

    matlist `diagnostics', title("SPUR diagnostics")
    matlist `levels', title("Levels SCPC")
    matlist `adjusted', title("Transformed SCPC")
    local row = 0
    foreach test in i0 i1 i0resid i1resid {
        local ++row
        local alldir = `info'[`row',2]
        local power = `info'[`row',3]
        _spatial_spur advisory, q(`q') nrep(`nrep') alldir(`alldir') power(`power') testname("SPUR `test'")
        return scalar `test'_p = `diagnostics'[`row',2]
        return scalar `test'_teststat = `diagnostics'[`row',1]
        return scalar `test'_ha_param = `diagnostics'[`row',3]
    }
    di as text "Note: LBM-GLS transformed variables are weighted-centered rather than untransformed."
    return matrix diagnostics = `diagnostics'
    return matrix diagnostic_info = `info'
    return matrix diagnostic_cutoffs = `cuts'
    return matrix levels_scpcstats = `levels'
    return matrix transformed_scpcstats = `adjusted'
    return matrix levels_critical_95 = `cv_levels'
    return matrix transformed_critical_95 = `cv_adjusted'
    if "`cvs'" != "" {
        return matrix levels_scpccvs = `extra_levels'
        return matrix transformed_scpccvs = `extra_adjusted'
    }
    return scalar q = `q'
    return scalar nrep = `nrep'
    return scalar tuning_draws = `nrep'
    return scalar n_obs = `n'
    return scalar n_sites = `n'
    return scalar levels_q = `levels_q'
    return scalar transformed_q = `transformed_q'
    if `"`drawfile'"' == "" return scalar seed = `seed'
    return local outputs "`outputs'"
    return local backend "rust"
    return local specification "`specification'"
    return local provider "`provider'"
end
