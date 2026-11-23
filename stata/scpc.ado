*! Spatial correlation-robust inference
program scpc, eclass sortpreserve
    version 16.0
    syntax , [CVS K(integer -1) AVC(real .03) LATlong UNCond EXACT APPROX METHOD(string) ///
        SEED(real 1) LANDmarks(integer 1000) SAMples(integer 20) TOLerance(real .001) ///
        CLuster(varname) BACKend(string) NThreads(integer -1) WEIGHTS CAPabilities]
    if "`capabilities'" != "" exit

    local estimator "`e(cmd)'"
    if !inlist("`estimator'", "regress", "areg", "ivregress") {
        di as error "SCPC supports unweighted regress, one-way areg, and ivregress 2sls"
        exit 198
    }
    if "`e(wtype)'" != "" {
        di as error "SCPC requires an unweighted fit"
        exit 198
    }
    if "`estimator'" == "ivregress" {
        if "`e(absvar)'" != "" {
            di as error "SCPC does not support ivregress with absorb()"
            exit 198
        }
        if "`e(estimator)'" != "2sls" {
            di as error "SCPC supports ordinary ivregress 2sls"
            exit 198
        }
    }
    if "`estimator'" == "areg" {
        local fe "`e(absvar)'"
        local count : word count `fe'
        if `count' != 1 {
            di as error "SCPC supports areg with one plain categorical absorb() variable"
            exit 198
        }
        capture confirm numeric variable `fe', exact
        if _rc {
            di as error "SCPC supports areg with one plain categorical absorb() variable"
            exit 198
        }
        if "`uncond'" == "" {
            di as error "SCPC after areg requires unconditional inference"
            exit 198
        }
    }
    if ("`exact'" != "" & "`approx'" != "") | ///
        ("`method'" != "" & ("`exact'" != "" | "`approx'" != "")) {
        di as error "Specify the geometry method once"
        exit 198
    }
    if "`exact'" != "" local method exact
    if inlist("`approx'", "approx") & "`approx'" != "" local method nystrom
    if "`method'" == "approx" local method nystrom
    if "`method'" == "" local method auto
    if !inlist("`method'", "auto", "exact", "nystrom") | `k' == 0 | `k' < -1 {
        di as error "Use method(auto), method(exact), or method(nystrom) and a positive k()"
        exit 198
    }
    if !inrange(`seed',0,4294967295) | `seed' != floor(`seed') {
        di as error "seed() must be an integer in [0, 2^32)"
        exit 198
    }

    _spatial_spur load, backend(`"`backend'"') nthreads(`nthreads')
    local provider "`r(provider)'"
    local specification "`r(specification)'"
    tempvar esample residual
    quietly gen byte `esample' = e(sample)
    _spatial_spur coords, sample(`esample') `latlong' minobs(5)
    local coords `r(coords)'
    local d = r(d)
    local n_sites = r(n)
    if "`cluster'" == "" {
        local cluster `e(clustvar)'
        if `: word count `cluster'' > 1 {
            di as error "Supply one actual spatial site identifier using cluster()"
            exit 198
        }
    }
    local groups
    if "`cluster'" != "" {
        tempvar group_id
        quietly count if `esample' & missing(`cluster')
        if r(N) {
            di as error "Cluster labels must be present for all fitted observations"
            exit 198
        }
        quietly egen long `group_id' = group(`cluster') if `esample'
        quietly summarize `group_id', meanonly
        local n_sites = r(max)
        local groups `group_id'
    }
    if "`method'" == "auto" local method = cond(`n_sites' >= 4500, "nystrom", "exact")
    quietly predict double `residual' if `esample', residuals

    tempname estimates omitted subset
    matrix `estimates' = e(b)
    quietly _ms_omit_info `estimates'
    matrix `omitted' = r(omit)
    local columns : colnames `estimates'
    local indices
    local actual_names
    local regressors
    local has_cons = 0
    local j = 0
    local endog_count = 0
    local endog `"`e(instd)'"'
    local endog : subinstr local endog "bn." ".", all
    local exog_terms
    foreach term of local columns {
        local ++j
        if `omitted'[1,`j'] != 0 continue
        if "`term'" == "_cons" {
            if "`estimator'" == "areg" continue
            local has_cons = 1
        }
        else {
            local normalized : subinstr local term "bn." ".", all
            local is_endog : list normalized in endog
            if `is_endog' {
                if "`exog_terms'" != "" {
                    di as error "Fitted endogenous and exogenous columns cannot be aligned"
                    exit 198
                }
                local ++endog_count
            }
            else local exog_terms `exog_terms' `normalized'
            quietly fvrevar `term' if `esample'
            local materialized `r(varlist)'
            if `: word count `materialized'' != 1 {
                di as error "Cannot identify fitted regressor `term'"
                exit 198
            }
            local regressors `regressors' `materialized'
        }
        local indices `indices' `j'
        local actual_names `actual_names' `term'
    }
    local p : word count `indices'
    if `p' == 0 {
        di as error "SCPC requires at least one estimable coefficient"
        exit 198
    }
    local report = cond(`k' == -1, `p', min(`k', `p'))
    mata: st_matrix("`subset'", st_matrix("`estimates'")[1, strtoreal(tokens("`indices'"))])

    local excluded
    if "`estimator'" == "ivregress" {
        if `endog_count' == 0 {
            di as error "Cannot identify fitted endogenous regressors"
            exit 198
        }
        local instruments `"`e(insts)'"'
        foreach term of local instruments {
            quietly _ms_parse_parts `term'
            if r(omit) continue
            local normalized : subinstr local term "bn." ".", all
            local is_exog : list normalized in exog_terms
            if `is_exog' continue
            quietly fvrevar `term' if `esample'
            local excluded `excluded' `r(varlist)'
        }
    }
    if "`estimator'" == "areg" {
        local within
        foreach var of local regressors {
            tempvar mean centered
            quietly egen double `mean' = mean(`var') if `esample', by(`fe')
            quietly gen double `centered' = `var' - `mean' if `esample'
            local within `within' `centered'
        }
        local regressors `within'
    }

    local nz : word count `excluded'
    local ll = "`latlong'" != ""
    local un = "`uncond'" != ""
    local cv = "`cvs'" != ""
    local cl = "`groups'" != ""
    local seed_arg = trim(string(`seed',"%15.0f"))
    local avc_arg = trim(string(`avc',"%21.16g"))
    local tolerance_arg = trim(string(`tolerance',"%21.16g"))
    tempname table cuts crit meta w site_rows
    matrix `table' = J(`report',6,.)
    matrix `crit' = J(`report',1,.)
    matrix `meta' = J(1,8,.)
    local cuts_arg _
    if `cv' {
        matrix `cuts' = J(`report',4,.)
        local cuts_arg `cuts'
    }
    tempfile weights_data
    local weights_path _
    if "`weights'" != "" local weights_path `"`weights_data'"'
    _spatial_spur invoke, vars("`coords' `residual' `groups' `regressors' `excluded'") sample(`esample') ///
        command(`"fitted_scpc `d' `p' `report' `has_cons' `endog_count' `nz' `cl' `ll' `avc_arg' `un' `cv' "`nthreads'" `method' `seed_arg' `landmarks' `samples' `tolerance_arg' `subset' `table' `cuts_arg' `crit' `meta' "`weights_path'""')

    local labels
    forvalues j = 1/`report' {
        local name : word `j' of `actual_names'
        local labels `labels' `name'
    }
    matrix rownames `table' = `labels'
    matrix colnames `table' = coef se t p ci_lower ci_upper
    matrix rownames `crit' = `labels'
    matrix colnames `crit' = critical_95
    if `cv' {
        matrix rownames `cuts' = `labels'
        matrix colnames `cuts' = cv32 cv10 cv05 cv01
    }
    if "`weights'" != "" {
        mata: spatial_spur_read_weights("`w'", "`site_rows'", `"`weights_path'"', ///
            st_matrix("`meta'")[1,7], st_matrix("`meta'")[1,4]+1, "`esample'", "`groups'")
        ereturn matrix W = `w'
        ereturn matrix site_rows = `site_rows'
    }
    if `cv' ereturn matrix scpccvs = `cuts'
    local inference = cond(`un', "unconditional", cond(`endog_count', "conditional_2sls", "conditional_ols"))
    matlist `table', title("SCPC `inference', `method'")
    ereturn matrix scpcstats = `table'
    ereturn matrix critical_95 = `crit'
    ereturn scalar c0 = `meta'[1,1]
    ereturn scalar cmax = `meta'[1,2]
    ereturn scalar cv = `meta'[1,3]
    ereturn scalar q = `meta'[1,4]
    ereturn scalar tied = `meta'[1,5]
    ereturn scalar n_obs = `meta'[1,6]
    ereturn scalar n_sites = `meta'[1,7]
    ereturn scalar distance_scale = `meta'[1,8]
    ereturn scalar avc = `avc'
    ereturn scalar cv_tolerance = `tolerance'
    ereturn scalar nthreads_requested = `nthreads'
    if "`method'" == "nystrom" {
        ereturn scalar seed = `seed'
        ereturn scalar landmarks = min(`landmarks', `n_sites')
        ereturn scalar landmark_samples = `samples'
    }
    local units = cond(`ll', "fraction of Earth circumference", "coordinate units")
    ereturn local inference "`inference'"
    ereturn local method "`method'"
    ereturn local backend "rust"
    ereturn local specification "`specification'"
    ereturn local provider "`provider'"
    ereturn local spatial_cluster "`cluster'"
    ereturn local coordinate_units "`units'"
end

mata:
void spatial_spur_read_weights(string scalar weights, string scalar row_matrix,
    string scalar path, real scalar n, real scalar q, string scalar sample, string scalar groups) {
    real rowvector representatives
    real colvector rows, ids
    real colvector format
    real scalar i, fd
    rows = selectindex(st_data(., sample))
    if (groups == "") representatives = rows'
    else {
        ids = st_data(rows, groups)
        representatives = J(1, n, .)
        for (i=1; i<=rows(rows); i++) {
            if (missing(representatives[ids[i]])) representatives[ids[i]] = rows[i]
        }
    }
    fd = fopen(path, "r")
    format = bufio()
    bufbyteorder(format, 2)
    st_matrix(weights, rowshape(fbufget(format, fd, "%8z", 1, n*q), n))
    fclose(fd)
    st_matrix(row_matrix, representatives)
}
end
