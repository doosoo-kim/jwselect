*! jwselect 0.9.0 04nov2026: scalar target adapter with full joint influence covariance
program define _jwselect_scalar, rclass
    version 16.1
    syntax [, VCE(string)]
    if !inlist("`vce'", "", "supplied") {
        di as err "Scalar target sources require vce(supplied) with their complete influence covariance"
        exit 198
    }
    if "`e(jw_schema)'" != "scalar_v1" {
        di as err "The source does not provide the validated scalar-target schema"
        exit 498
    }
    if `"`e(datasignaturevars)'"' == "" {
        di as err "Scalar target sources must provide a signed estimation sample; refit its source"
        exit 498
    }
    capture checkestimationsample
    if _rc {
        di as err "Data changed since the scalar target was fitted; refit its source"
        exit 459
    }
    tempname d V bases target h hm
    matrix `d' = e(jw_d)
    matrix `V' = e(jw_V)
    matrix `bases' = e(jw_bases)
    matrix `target' = e(jw_targetmap)
    local K = e(jw_K)
    if missing(`K') | `K'<0 | `K'!=floor(`K') error 498
    if `K' {
        matrix `h' = e(jw_h)
        matrix `hm' = e(jw_hmap)
        local delta=e(anticipation)
        local gap=e(gap)
        if missing(`delta') | `delta'<0 | `delta'!=floor(`delta') error 498
        forvalues j=1/`K' {
            if `hm'[`j',2] > `hm'[`j',1]-(`delta'+1)*`gap' {
                di as err "An optional restriction is inside the source anticipation window"
                exit 498
            }
        }
        return matrix h = `h'
        return matrix hmap = `hm'
    }
    return matrix d = `d'
    return matrix V = `V'
    return matrix bases = `bases'
    return matrix targetmap = `target'
    return scalar K = `K'
    return scalar N_units = e(jw_nunits)
    return scalar N_obs = e(N)
    return scalar N_clusters = e(jw_nunits)
    return scalar N_covariates = e(N_covariates)
    return scalar gap = e(gap)
    return scalar anticipation = e(anticipation)
    return local target_label `"`e(jw_target_label)'"'
    return local restriction_scale `"`e(jw_restriction_scale)'"'
    return local sourcecmd `"`e(cmdline)'"'
    return local method `"`e(method)'"'
    return local model scalar_target
    return local sourceversion 0.9.0
    return local upstream_stat z
end
