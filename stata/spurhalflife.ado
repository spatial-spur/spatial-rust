*! Spatial half-life confidence set
program spurhalflife, rclass sortpreserve
    version 16.0
    syntax varname(numeric) [if] [in], [Q(integer 15) NRep(integer 100000) SEED(real 42) ///
        LEVEL(real 95) LATlong NORMdist BACKend(string) NThreads(integer -1) DRAWfile(string)]
    if `level' <= 0 | `level' >= 100 | `q' < 2 | `q' > 100 | `nrep' < 1 {
        di as error "Half-life requires 2 <= q <= 100, positive nrep(), and 0 < level < 100"
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
    if `q' >= `n' {
        di as error "q() must be smaller than the sample"
        exit 198
    }
    local use_file `"`drawfile'"'
    if `"`use_file'"' == "" {
        tempfile normals
        _spatial_spur draws, q(`q') nrep(`nrep') seed(`seed') file(`"`normals'"')
        local use_file `"`normals'"'
    }
    tempname lower upper diameter info components count
    matrix `info' = J(1,3,.)
    matrix `components' = J(65,2,.)
    local ll = "`latlong'" != ""
    _spatial_spur invoke, vars("`coords' `varlist'") sample(`touse') ///
        command(`"spur_halflife `d' `ll' `q' `nrep' `level' "`nthreads'" "`use_file'" `lower' `upper' `diameter' `info' `components' `count'"')
    local count_value = scalar(`count')
    local normal = "`normdist'" != ""
    if !`normal' {
        scalar `lower' = scalar(`lower') * scalar(`diameter')
        scalar `upper' = scalar(`upper') * scalar(`diameter')
        mata: st_matrix("`components'", st_matrix("`components'") * st_numscalar("`diameter'"))
    }
    local units = cond(`ll', "meters", "coordinate units")
    local display_units = cond(`normal', "sample diameters", "`units'")
    di as text "`level'% spatial half-life confidence set (`display_units')"
    if `count_value' > 0 {
        matrix `components' = `components'[1..`count_value',1..2]
        matrix colnames `components' = lower upper
        forvalues i = 1/`count_value' {
            local right = cond(missing(`components'[`i',2]), "Inf", trim(string(`components'[`i',2], "%12.6g")))
            di as result "[" %12.6g `components'[`i',1] ", `right']"
        }
        return matrix components = `components'
    }
    else di as result "Empty confidence set"
    local alldir = `info'[1,2]
    _spatial_spur advisory, q(`q') nrep(`nrep') level(`level') alldir(`alldir')
    return scalar ci_l = scalar(`lower')
    return scalar ci_u = scalar(`upper')
    return scalar max_dist = scalar(`diameter')
    return scalar n_components = `count_value'
    return scalar empty = `count_value' == 0
    return scalar unbounded = `count_value' > 0 & missing(scalar(`upper'))
    return scalar tied = `info'[1,1]
    return scalar all_directions = `info'[1,2]
    return scalar n_obs = `n'
    return scalar n_sites = `n'
    return scalar q = `q'
    return scalar nrep = `nrep'
    return scalar level = `level'
    return scalar normdist = `normal'
    if `"`drawfile'"' == "" return scalar seed = `seed'
    return local distance_units "`units'"
    return local normalization "sample diameter"
    return local backend "rust"
    return local specification "`specification'"
    return local provider "`provider'"
end
