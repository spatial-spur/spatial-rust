*! Spatial SPUR shared Stata adapter
program _spatial_spur, rclass
    version 16.0
    gettoken operation 0 : 0, parse(" ,")

    if "`operation'" == "load" {
        syntax [, BACKend(string) nthreads(integer -1)]
        local backend = lower(trim(`"`backend'"'))
        if !inlist("`backend'", "", "auto", "rust") {
            di as error "backend() must be rust or auto"
            exit 198
        }
        if `nthreads' < 1 & `nthreads' != -1 {
            di as error "nthreads() must be positive"
            exit 198
        }
        capture plugin call _spatial_spur_plugin, ping_unified_v1
        if _rc {
            capture findfile spatial_spur.plugin
            if _rc {
                di as error "Install the matching spatial_spur.plugin for this Stata platform"
                exit 199
            }
            capture program _spatial_spur_plugin, plugin using("`r(fn)'")
            capture plugin call _spatial_spur_plugin, ping_unified_v1
            if _rc {
                di as error "The installed spatial_spur.plugin is incompatible"
                exit 199
            }
        }
        capture plugin call _spatial_spur_plugin, require_spec 0.5.0
        if _rc {
            di as error "The installed spatial_spur.plugin does not match specification 0.5.0"
            exit 199
        }
        quietly plugin call _spatial_spur_plugin, result_metadata
        return local backend "rust"
        return local provider "`native_provider'"
        return local specification "`native_methods'"
        exit
    }

    if "`operation'" == "invoke" {
        syntax , VARS(string) SAMPLE(varname) COMMAND(string)
        plugin call _spatial_spur_plugin `vars' if `sample', `command'
        exit
    }

    if "`operation'" == "coords" {
        syntax , SAMPLE(varname) [latlong MINobs(integer 2)]
        quietly count if `sample'
        local n = r(N)
        if `n' < `minobs' {
            di as error "Too few observations in the selected sample"
            exit 2001
        }
        capture unab vars : s_*
        if _rc {
            di as error "Supply coordinates in s_1, s_2, and so on"
            exit 111
        }
        local d : word count `vars'
        local coords
        forvalues j = 1/`d' {
            capture confirm numeric variable s_`j', exact
            if _rc {
                di as error "Spatial coordinates must be consecutively named numeric variables starting with s_1"
                exit 198
            }
            quietly count if `sample' & missing(s_`j')
            if r(N) > 0 {
                di as error "Coordinate s_`j' contains missing fitted observations"
                exit 198
            }
            local coords `coords' s_`j'
        }
        if "`latlong'" != "" {
            if `d' != 2 {
                di as error "latlong requires s_1 latitude and s_2 longitude"
                exit 198
            }
            quietly count if `sample' & !inrange(s_1,-90,90)
            if r(N) > 0 {
                di as error "Latitude s_1 must be between -90 and 90 degrees"
                exit 198
            }
        }
        return scalar n = `n'
        return scalar d = `d'
        return local coords "`coords'"
        exit
    }

    if "`operation'" == "draws" {
        syntax , Q(integer) NRep(integer) SEED(real) FILE(string)
        if `q' < 1 | `nrep' < 1 | !inrange(`seed',0,2147483647) | `seed' != floor(`seed') {
            di as error "Draw dimensions must be positive and seed() must be an integer in [0, 2147483647]"
            exit 198
        }
        local seed_str = trim(string(`seed', "%15.0f"))
        local rngstate `"`c(rngstate)'"'
        capture erase `"`file'"'
        tempname m_fd m_fmt m_draws
        mata: `m_fd' = -1
        capture noisily {
            quietly set seed `seed_str'
            mata: `m_draws' = rnormal(`q', `nrep', 0, 1); `m_fmt' = bufio(); bufbyteorder(`m_fmt', 2); `m_fd' = fopen(st_local("file"), "w"); fbufput(`m_fmt', `m_fd', "%8z", rowshape(`m_draws', 1)); fclose(`m_fd'); `m_fd' = -1; `m_draws' = .
        }
        local rc = _rc
        capture mata: if (`m_fd' >= 0) fclose(`m_fd')
        capture mata: mata drop `m_fd' `m_fmt' `m_draws'
        quietly set rngstate `rngstate'
        if `rc' exit `rc'
        exit
    }

    if "`operation'" == "advisory" {
        syntax , Q(integer) NRep(integer) [LEVEL(real -1) ALLdir(integer 0) POWer(real -1) TESTname(string)]
        if `q' >= 40 | `alldir' == 1 {
            local scope = cond(`alldir' == 1, " uses all positive directions and", "")
            di as text "Note: q=`q'`scope' retains higher spatial frequencies. Results can be sensitive to local noise."
        }
        if `level' < 0 {
            local unatt
            foreach alpha in 1 5 10 {
                if `nrep' + 1 < floor(100 / `alpha') {
                    local unatt `unatt' `alpha'%
                }
            }
            if "`unatt'" != "" {
                local minp = trim(string(1 / (`nrep' + 1), "%8.4g"))
                di as text "Note: With nrep=`nrep', minimum p-value is `minp'. Rejection is unattainable at`unatt'."
            }
            if `power' >= 0 & `power' < 0.45 {
                local pct = trim(string(100 * `power', "%5.1f"))
                di as text "Note: `testname' did not reach the 50% tuning target (best estimated power `pct'%)."
            }
        }
        else {
            if `q' <= 5 {
                di as text "Note: q=`q' provides few spatial contrasts for estimating the half-life."
            }
            if (`nrep' / (`nrep' + 1)) < (`level' / 100) {
                local minp = trim(string(1 / (`nrep' + 1), "%8.4g"))
                di as text "Note: With nrep=`nrep', minimum p-value is `minp'. The `level'% confidence set covers [0, infinity]."
            }
        }
        exit
    }

    if inlist("`operation'", "validate", "commit") {
        syntax , OUTPUTS(string) [SOURCE(string) REPLACE]
        local seen
        foreach name of local outputs {
            if `: list name in seen' {
                di as error "Duplicate output variable name `name'"
                exit 198
            }
            local seen `seen' `name'
            capture confirm name `name'
            if _rc | strlen("`name'") > 32 {
                di as error "Invalid output variable name `name'"
                exit 198
            }
            capture confirm variable `name', exact
            if !_rc & "`replace'" == "" {
                di as error "Output variable `name' already exists"
                exit 110
            }
        }
        if "`operation'" == "commit" {
            preserve
            foreach name of local outputs {
                capture confirm variable `name', exact
                if !_rc quietly drop `name'
            }
            rename (`source') (`outputs')
            restore, not
        }
        exit
    }

    di as error "Unknown spatialSpur adapter operation `operation'"
    exit 198
end
