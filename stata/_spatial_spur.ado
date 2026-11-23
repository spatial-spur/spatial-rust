*! Spatial SPUR shared Stata adapter
program _spatial_spur, rclass
    version 16.0
    gettoken operation 0 : 0, parse(" ,")

    if "`operation'" == "load" {
        syntax [, BACKend(string) nthreads(integer -1)]
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
        quietly plugin call _spatial_spur_plugin, result_metadata
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
        if `q' < 1 | `nrep' < 1 | !inrange(`seed',0,4294967295) | `seed' != floor(`seed') {
            di as error "Draw dimensions and seed() must be valid integer values"
            exit 198
        }
        local rngstate `"`c(rngstate)'"'
        capture noisily {
            quietly set seed `seed'
            mata: spatial_spur_write_draws(`"`file'"', `q', `nrep')
        }
        local rc = _rc
        quietly set rngstate `rngstate'
        if `rc' exit `rc'
        exit
    }

    if inlist("`operation'", "validate", "commit") {
        syntax , OUTPUTS(string) [SOURCE(string) REPLACE]
        foreach name of local outputs {
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

mata:
void spatial_spur_write_draws(string scalar path, real scalar q, real scalar nrep) {
    real matrix draws
    real scalar fd
    real colvector format
    draws = rnormal(q, nrep, 0, 1)
    fd = fopen(path, "w")
    format = bufio()
    bufbyteorder(format, 2)
    fbufput(format, fd, "%8z", rowshape(draws, 1))
    fclose(fd)
}
end
