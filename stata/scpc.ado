*! Spatial correlation-robust inference
program scpc, eclass sortpreserve
    version 16.0
    syntax , [CVS K(integer -1) AVC(real .03) LATlong UNCond EXACT APPROX METHOD(string) ///
        SEED(real 1) LANDmarks(integer 1000) SAMples(integer 20) TOLerance(real .001) ///
        CLuster(varname) BACKend(string) NThreads(integer -1) WEIGHTS CAPabilities]
    if "`capabilities'" != "" {
        _spatial_spur load, backend(`"`backend'"') nthreads(`nthreads')
        di as text "specification: " as result "`r(specification)'" ///
            as text "  backend: " as result "`r(backend)'" ///
            as text "  provider: " as result "`r(provider)'"
        exit
    }

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
        capture confirm variable `fe', exact
        if _rc {
            di as error "SCPC supports areg with one plain categorical absorb() variable"
            exit 198
        }
        if "`uncond'" == "" {
            di as error "SCPC after areg requires unconditional inference"
            exit 198
        }
    }
    local method = lower(trim(`"`method'"'))
    if ("`exact'" != "" & "`approx'" != "") | ///
        ("`method'" != "" & ("`exact'" != "" | "`approx'" != "")) {
        di as error "Specify the geometry method once"
        exit 198
    }
    if "`exact'" != "" local method exact
    if "`approx'" != "" local method nystrom
    if "`method'" == "approx" local method nystrom
    if "`method'" == "" local method auto
    if !inlist("`method'", "auto", "exact", "nystrom") | `k' == 0 | `k' < -1 {
        di as error "Use method(auto), method(exact), or method(nystrom) and a positive k()"
        exit 198
    }
    if !(`avc' > 0.001 & `avc' < 0.99) {
        di as error "avc() must be strictly between 0.001 and 0.99"
        exit 198
    }
    if !(`tolerance' > 0 & `tolerance' < .) | `landmarks' < 2 | `samples' < 1 {
        di as error "Require tolerance() > 0, landmarks() >= 2, and samples() >= 1"
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
    tempvar squared_residual
    quietly gen double `squared_residual' = `residual'^2 if `esample'
    quietly summarize `squared_residual', meanonly
    if r(N) != e(N) | reldif(r(sum), e(rss)) > 1e-9 {
        di as error "The current observations no longer reproduce the fitted regression residuals"
        exit 198
    }


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
        if `has_cons' {
            di as error "_cons must be the last fitted column"
            exit 198
        }
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
            local materialized `r(varlist)'
            if `: word count `materialized'' != 1 {
                di as error "Cannot identify excluded instrument `term'"
                exit 198
            }
            local excluded `excluded' `materialized'
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
    mata: st_matrix("e(scpccvs)", J(0, 0, .))
    mata: st_matrix("e(W)", J(0, 0, .))
    mata: st_matrix("e(site_rows)", J(0, 0, .))
    mata: st_numscalar("e(seed)", J(0, 0, .))
    mata: st_numscalar("e(landmarks)", J(0, 0, .))
    mata: st_numscalar("e(landmark_samples)", J(0, 0, .))

    if "`weights'" != "" {
        local w_n = `meta'[1,7]
        local w_q = `meta'[1,4] + 1
        if `w_n' > c(max_matsize) | `w_q' > c(max_matsize) {
            di as error "weights matrix dimensions (`w_n' x `w_q') exceed max_matsize (`c(max_matsize)')"
            exit 198
        }
        tempname m_fd m_fmt m_rows m_ids m_rep m_i
        mata: `m_fd' = -1
        capture noisily {
            mata: `m_rows' = selectindex(st_data(., "`esample'"))
            if "`groups'" == "" {
                mata: `m_rep' = `m_rows''
            }
            else {
                mata: `m_ids' = st_data(`m_rows', "`groups'")
                mata: `m_rep' = J(1, `w_n', .)
                mata: for (`m_i' = 1; `m_i' <= rows(`m_rows'); `m_i'++) if (missing(`m_rep'[`m_ids'[`m_i']])) `m_rep'[`m_ids'[`m_i']] = `m_rows'[`m_i']
            }
            mata: `m_fd' = fopen(st_local("weights_path"), "r"); `m_fmt' = bufio(); bufbyteorder(`m_fmt', 2); st_matrix("`w'", rowshape(fbufget(`m_fmt', `m_fd', "%8z", 1, `w_n' * `w_q'), `w_n')); fclose(`m_fd'); `m_fd' = -1
            mata: st_matrix("`site_rows'", `m_rep')
        }
        local rc = _rc
        capture mata: if (`m_fd' >= 0) fclose(`m_fd')
        capture mata: mata drop `m_fd' `m_fmt' `m_rows' `m_ids' `m_rep' `m_i'
        if `rc' exit `rc'
        ereturn matrix W = `w'
        ereturn matrix site_rows = `site_rows'
    }
    local inference = cond(`un', "unconditional", cond(`endog_count', "conditional_2sls", "conditional_ols"))
    matlist `table', title("SCPC `inference', `method'")
    local summary_n_obs = `meta'[1,6]
    local summary_n_sites = `meta'[1,7]
    local summary_q = `meta'[1,4]
    di as text "n_obs = " as result `summary_n_obs' ///
        as text "  n_sites = " as result `summary_n_sites' ///
        as text "  q = " as result `summary_q' ///
        as text "  avc = " as result %6.4g `avc'
    if "`method'" == "nystrom" {
        local used_landmarks = min(`landmarks', `n_sites')
        di as text "Nystrom: seed = " as result `seed' ///
            as text "  landmarks = " as result `used_landmarks' ///
            as text "  samples = " as result `samples'
    }
    if `meta'[1,5] == 1 {
        di as text "Note: Spatial correlation matrix had tied boundary eigenvalues; critical values reflect a randomized eigenspace basis."
    }
    if "`inference'" == "conditional_2sls" {
        di as text "Note: Conditional 2SLS SCPC inference assumes strong first-stage identification; under weak instruments use Anderson-Rubin inference on reduced-form projections."
    }
    matlist `crit', title("SCPC 95% critical values")
    if `cv' {
        matlist `cuts', title("SCPC critical values")
        ereturn matrix scpccvs = `cuts'
    }
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
