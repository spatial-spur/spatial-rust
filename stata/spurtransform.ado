*! Spatial transformation with staged output variables
program spurtransform, rclass sortpreserve
    version 16.0
    syntax varlist(numeric) [if] [in], PREFIX(name) [TRANSformation(string) RADIUS(string) ///
        CLUSTVAR(varname) LATlong REPLACE SEParately BACKend(string) NThreads(integer -1)]
    if "`transformation'" == "" local transformation lbmgls
    if !inlist("`transformation'", "lbmgls", "nn", "iso", "cluster") {
        di as error "Choose lbmgls, nn, iso, or cluster"
        exit 198
    }
    if ("`transformation'" == "iso") != (trim("`radius'") != "") {
        di as error "Supply radius() only for ISO"
        exit 198
    }
    if ("`transformation'" == "cluster") != ("`clustvar'" != "") {
        di as error "Supply clustvar() only for the cluster transform"
        exit 198
    }
    local outputs
    foreach v of local varlist {
        local outputs `outputs' `prefix'`v'
    }
    _spatial_spur validate, outputs("`outputs'") `replace'
    _spatial_spur load, backend(`"`backend'"') nthreads(`nthreads')
    local provider "`r(provider)'"
    local specification "`r(specification)'"

    local source
    foreach v of local varlist {
        tempvar result
        quietly gen double `result' = .
        local source `source' `result'
    }
    local coords
    local groups
    local d 0
    local ll = "`latlong'" != ""
    local cluster = "`clustvar'" != ""
    if `cluster' {
        tempvar group
        quietly egen long `group' = group(`clustvar')
        local groups `group'
    }
    if "`radius'" == "" local radius -1
    if "`separately'" == "" {
        marksample touse
        if `cluster' quietly replace `touse' = 0 if missing(`groups')
        else {
            _spatial_spur coords, sample(`touse') `latlong'
            local coords `r(coords)'
            local d = r(d)
        }
        local variables : word count `varlist'
        _spatial_spur invoke, vars("`coords' `groups' `varlist' `source'") sample(`touse') ///
            command(`"spur_transform `d' `transformation' `ll' "`radius'" `cluster' `variables' "`nthreads'""')
    }
    else {
        marksample selected, novarlist
        local i = 0
        foreach v of local varlist {
            local ++i
            tempvar touse
            quietly gen byte `touse' = `selected' & !missing(`v')
            if `cluster' quietly replace `touse' = 0 if missing(`groups')
            else {
                _spatial_spur coords, sample(`touse') `latlong'
                local coords `r(coords)'
                local d = r(d)
            }
            local result : word `i' of `source'
            _spatial_spur invoke, vars("`coords' `groups' `v' `result'") sample(`touse') ///
                command(`"spur_transform `d' `transformation' `ll' "`radius'" `cluster' 1 "`nthreads'""')
        }
    }
    _spatial_spur commit, outputs("`outputs'") source("`source'") `replace'
    return local outputs "`outputs'"
    return local transformation "`transformation'"
    return local backend "rust"
    return local specification "`specification'"
    return local provider "`provider'"
end
